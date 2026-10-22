#!/usr/bin/env python3
"""
Reorganize the azcot hourly climatology archives into a single
timeseries-optimized Zarr store.

Source can be a local filesystem path (e.g. a Qumulo NFS/SMB mount) or
an HTTP root - auto-detected from the --root value. Default is the
local mount:

    /qumulo/azcot/2T/<MM>/<DD>/2T.<mon><DD><HH>.91-20clim.nc
    /qumulo/azcot/rad/<MM>/<MM><DD><HH>.91-20radclim.nc
    /qumulo/azcot/windclimo/<MM>/<DD>/<MM><DD><HH>.WS10.nc

Reading local files skips the HTTP-specific workarounds entirely:
xr.open_dataset(path) opens a local file natively and fast (no need
to dodge OPeNDAP auto-detection or route bytes through fsspec, both of
which only matter for "http://" URLs - see open_source() below).

Why rechunk (not just repackage)
--------------------------------
Each source archive is one netCDF file per hour holding a full global
grid. That is optimal for "global map at one instant" and pessimal for
"timeseries at one point": a point-year needs 8,784 file reads and
several GB read to extract 8,784 floats - true whether those reads are
local disk I/O or HTTP requests.

Merging into one big netCDF/HDF5 file would not fix this for HTTP
consumers - without OPeNDAP the client must fetch the whole file. Zarr
splits the array into many small independently-addressable objects, so
a client (local or remote) reads only the chunks it needs.

The three source layouts are all different
------------------------------------------
Confirmed from the archive listings:

  2T    2T/<MM>/<DD>/2T.<mon><DD><HH>.91-20clim.nc     688K  day dirs, month ABBREV
  rad   rad/<MM>/<MM><DD><HH>.91-20radclim.nc          4.2M  FLAT, no day dirs
  wind  windclimo/<MM>/<DD>/<MM><DD><HH>.WS10.nc       689K  day dirs, numeric month

Notes on each:
  * 2T and WS10 are ~697 KB = one 121x1440 float32 field (0.25 deg).
  * WS10 is wind SPEED already computed - no u/v combination needed.
    (WD10, wind direction, sits alongside it and is ignored here. If you
    ever want direction, add it as another VarSpec, but note direction
    must NOT be averaged or interpolated linearly across the 0/360
    wrap - treat it separately.)
  * rad is 4.2 MB = one 301x3600 float32 field (0.10 deg) - i.e. a
    FINER GRID, not several variables. Same domain (lat 90..60,
    lon -180..180), 2.5x the resolution.
  * Confirmed on the real archive: rad's variable is literally named
    "SSRD_GDS0_SFC_acc1h" with units "W m**-2 s". Dimensionally
    W*s = J, so this is accumulated ENERGY (J m-2) over the named
    period, not instantaneous irradiance, despite "W" in the units
    string. See infer_deaccum_seconds() - the conversion to W m-2 is
    applied automatically (--no-deaccumulate to disable).

Two grids, two groups - no regridding
-------------------------------------
2T/WS10 (0.25 deg) and SSRD (0.10 deg) do not share a mesh, and this
script does NOT interpolate them onto one:

  * 0.25/0.10 = 2.5 is not an integer ratio, so cells coincide only
    every 0.5 deg. Any merge is genuine interpolation.
  * Downsampling SSRD throws away resolution the provider paid to
    produce; upsampling 2T/WS10 inflates them 6.25x for no new
    information.
  * A flux like solar irradiance wants area-weighted/conservative
    remapping, which is the data owner's scientific call.
  * For point timeseries it is simply unnecessary - take the nearest
    cell in each grid.

So each grid becomes its own group in one store:

    grid_0p25/   time, lat(121), lon(1440), 2T, WS10
    grid_0p10/   time, lat(301), lon(3600), SSRD

Both groups carry an identical time axis, so they stay aligned.
Read them with:

    a = xr.open_zarr(url, group="grid_0p25")
    b = xr.open_zarr(url, group="grid_0p10")

Storage, full 8,784-hour year (uncompressed / ~zstd):
    2T    121x1440    6.1 GB / ~3.4 GB
    WS10  121x1440    6.1 GB / ~3.4 GB
    SSRD  301x3600   38.1 GB / ~21.2 GB   <- the finer grid dominates

Chunk geometry
--------------
For the full 8,784-hour year on a 121x1440 grid:

    chunk (t,lat,lon)   chunk size   n chunks   1-pt-year: reqs / bytes
    ------------------------------------------------------------------
    (1, 121, 1440)        697 KB       8,784      8,784  /  6.1 GB  <- today
    (8784, 8, 8)          2.25 MB      2,880          1  /  2.2 MB
    (2196, 8, 8)          563 KB      11,520          4  /  2.2 MB  <- default
    (744, 16, 16)         762 KB       8,640         12  /  9.1 MB

Default (2196, 8, 8): ~563 KB chunks, 11,520 objects per variable, and
a point-year costs 4 requests / 2.2 MB instead of 8,784 / 6.1 GB -
roughly 2,800x fewer bytes.

Tradeoff: this makes single-timestep global maps worse (a map touches
2,880 chunks). Keep the original per-hour netCDFs alongside, or publish
a second Zarr chunked (1, 121, 1440), if map access matters.

Usage
-----
    # 1. See what variables actually live in the rad files:
    python reorganize_to_zarr.py --list-vars

    # 2. Build the store (all three variables), reading local files:
    python reorganize_to_zarr.py --out /srv/data/azcot/azcot.zarr

    # 3. Explicit source root (local path or http(s) URL):
    python reorganize_to_zarr.py --out out.zarr --root /qumulo/azcot
    python reorganize_to_zarr.py --out out.zarr \\
        --root https://scil-data.sdsc.edu/data/azcot

    # 4. Subset, with an explicit rad variable name:
    python reorganize_to_zarr.py --out out.zarr --vars 2T WS10 \\
        --rad-var SSRD_GDS0_SFC

Adding variables to a store that already exists
-----------------------------------------------
Never re-run a build with the same --out to add variables: a build
creates the output from scratch and overwrites it. Use --append, which
leaves everything already in the store untouched.

The newer families (STL1, SD, SKT, WCT) arrive as tarballs. Extract each
family into its own directory under the source root; the members are
already laid out as <MM>/<file>, so no further path handling is needed:

    cd /qumulo/azcot
    for fam in STL1 sd skt wct; do
      mkdir -p "$fam"
      for t in archive/*"$fam".tar.gz archive/*"$fam".tgz; do
        [ -e "$t" ] && tar -xzf "$t" -C "$fam"
      done
    done

(GNU tar keeps everything inside the target directory: it strips a
leading "/" from absolute member names and skips members containing
"..", exiting non-zero if it skipped any - so a tarball with such names
fails loudly instead of writing elsewhere. Only the climatology
families need extracting: the 2t, wind and rad tarballs are
mostly raw per-year files, a different product that does not belong on
this placeholder-year time axis, and their climatology is already in
the store.)

    # what would go in, and which months each source lacks:
    python reorganize_to_zarr.py --list-vars --vars STL1 SD SKT WCT

    # add them:
    python reorganize_to_zarr.py --out /qumulo/azcot/azcot.zarr --append \\
        --vars STL1 SD SKT WCT --max-workers 4

    # record units once established (writes no data):
    python reorganize_to_zarr.py --out /qumulo/azcot/azcot.zarr --annotate \\
        --meta WCT.units=degF

HDF5-format source files (NetCDF-4) are read under a process-wide lock,
because opening them from several threads at once through the netCDF4
engine corrupts memory (a segfault or "free(): invalid size" partway
through a long run). Classic NetCDF-3 files are unaffected and stay
parallel, so --max-workers speeds up those but not HDF5 ones. To see
which kind a file is:

    python3 -c "print(open('/qumulo/azcot/wct/01/010100.91-20wctclim.nc','rb').read(8))"
    # b'\\x89HDF\\r\\n\\x1a\\n' = HDF5 / NetCDF-4     b'CDF\\x01' or b'CDF\\x02' = classic

If a run is interrupted (crash, Ctrl-C, lost mount):

  * Adding NEW variables: the store is unaffected. They are registered
    only once all have been written, so what was written is an invisible
    leftover. Re-run the same command: leftovers of the variables being
    added are deleted and rewritten. Add --resume to keep the ones that
    finished completely and redo only the unfinished one.
  * --replace: the variable is un-registered BEFORE its data is deleted,
    so an interruption leaves it absent from the store (never listed but
    half destroyed, which an earlier version allowed). Re-run to restore
    it. It is unavailable to readers while it is being rewritten.
  * --fill: nothing is ever deleted, so the array is valid at every
    moment; re-run to finish.

Months that arrive after a variable was added
---------------------------------------------
Extract the new tarball(s) into the family's directory, then --fill. It
writes only the months whose source files now exist and that the array's
months_with_data does not list, in place, without deleting anything:

    tar -xzf archive/novSTL1.tar.gz -C /qumulo/azcot/STL1
    python reorganize_to_zarr.py --out /qumulo/azcot/azcot.zarr --fill \\
        --vars STL1 SD SKT --max-workers 4

Run it BEFORE verify_zarr.py: until then the store holds NaN where source
files now exist, and the verifier will (correctly) report that as a
mismatch. Each month is checked against what the array already holds
first: the files' units must equal the recorded source_units and the grid
must equal the group's, or the month is refused and nothing is written. A
month that already has data is refused too (use --replace to refresh).
Arrays built before months_with_data existed (the original 2T, WS10 and
SSRD) need --months to say which months are empty.

Prefer --fill to --replace for this: it reads only the new months, and
cannot leave a variable missing.

Append refuses, and changes nothing, if: a variable is already present
(use --replace NAME to overwrite deliberately); the source has a month
the store's time axis lacks; or the store is not one this tool made. A
month the store has but the source lacks is left NaN and recorded in the
array's months_with_data attribute. A variable on a grid no existing
group matches gets a new group; one with the same spacing but a
different shape gets a disambiguated name (e.g. grid_0p25_121x720)
rather than reusing or overwriting the existing group.

Clients read the finished store however they like (local path or, if
you serve it, a URL) and need none of the per-file scaffolding:

    ds = xr.open_zarr("/srv/data/azcot/azcot.zarr", group="grid_0p25")
    ds["WS10"].sel(lat=38.9, lon=-77.0, method="nearest").plot()

Requires: xarray, zarr>=2.18, numcodecs, numpy, pandas, scipy, h5netcdf,
          h5py, and (only if --root is an http(s) URL) fsspec, aiohttp

    pip install -U 'zarr>=2.18' numcodecs xarray numpy pandas scipy \\
        h5netcdf h5py fsspec aiohttp

Note on zarr versions: zarr < 2.18 uses the numpy aliases np.PINF /
np.NINF, removed in NumPy 2.0, so an old zarr on a modern NumPy fails
with "AttributeError: `np.PINF` was removed in the NumPy 2.0 release."
This script works with zarr 2.18+ or zarr 3.x (see the compat shim
below) and enforces the 2.18 floor at import.
"""

import argparse
import collections
import concurrent.futures as cf
import contextlib
import datetime as dt
import os
import re
import threading

import numpy as np
import pandas as pd
import xarray as xr
import zarr

# Works with both zarr-python 2.x and 3.x via the small shim below.
#
# The floor is zarr 2.18: earlier 2.x releases use the numpy aliases
# np.PINF / np.NINF, which NumPy 2.0 removed, so zarr < 2.18 on
# NumPy >= 2 dies with
#     AttributeError: `np.PINF` was removed in the NumPy 2.0 release.
# zarr 2.18 is the final 2.x series and is NumPy-2 clean.
_ZV = tuple(int(x) for x in str(zarr.__version__).split(".")[:2])
ZARR2 = _ZV[0] == 2
if _ZV < (2, 18):
    raise ImportError(
        f"zarr {zarr.__version__} is too old. Use zarr >= 2.18 (final 2.x)\n"
        f"or zarr >= 3:\n"
        f"    pip install -U 'zarr>=2.18,<3'      # stay on 2.x\n"
        f"    pip install -U 'zarr>=3'            # or move to 3.x\n"
        f"zarr < 2.18 crashes with 'np.PINF was removed in the NumPy 2.0 "
        f"release' when used with NumPy >= 2.")


def _make_compressor(clevel, zarr_format):
    """Blosc/zstd compressor in whichever form this zarr version wants."""
    import numcodecs
    if ZARR2 or zarr_format == 2:
        return numcodecs.Blosc("zstd", clevel=clevel,
                               shuffle=numcodecs.Blosc.SHUFFLE)
    return zarr.codecs.BloscCodec(cname="zstd", clevel=clevel,
                                  shuffle=zarr.codecs.BloscShuffle.shuffle)


def _open_group(path, zarr_format):
    """zarr 2 has no zarr_format kwarg; it only ever writes v2."""
    if ZARR2:
        return zarr.open_group(path, mode="w")
    return zarr.open_group(path, mode="w", zarr_format=zarr_format)


def _make_array(grp, name, shape, chunks, dtype, compressor=None,
                fill_value=None):
    """
    zarr 2: Group.create_dataset(..., compressor=<numcodecs codec>)
    zarr 3: Group.create_array(...,  compressors=[<codec>])
    """
    if ZARR2:
        return grp.create_dataset(name, shape=shape, chunks=chunks,
                                  dtype=dtype, compressor=compressor,
                                  fill_value=fill_value)
    kw = {} if compressor is None else {"compressors": [compressor]}
    return grp.create_array(name, shape=shape, chunks=chunks, dtype=dtype,
                            fill_value=fill_value, **kw)

DEFAULT_ROOT = "/qumulo/azcot"   # local mount; pass --root for an http(s) URL

MONTH_ABBR = ["jan", "feb", "mar", "apr", "may", "jun",
              "jul", "aug", "sep", "oct", "nov", "dec"]

PLACEHOLDER_YEAR = 2000          # climatology has no real year; leap-safe
DEFAULT_CHUNKS = (2196, 8, 8)    # (time, lat, lon)


def is_url(root):
    return root.startswith("http://") or root.startswith("https://")

# Candidate names for surface solar irradiance inside the rad files,
# most-likely first. NOT verified against the actual files - the data
# host is unreachable from the machine this was written on.
# --list-vars will show the truth.
RAD_VAR_PREFERENCE = [
    "SSRD_GDS0_SFC",   # surface solar radiation downwards (ECMWF/NCL style)
    "SSRD",
    "ssrd",
    "SSR_GDS0_SFC",    # net surface solar (downwards minus reflected)
    "SSR",
]

# Matches an ECMWF-style accumulation-period tag on the variable name,
# e.g. "SSRD_GDS0_SFC_acc1h" -> 1. Confirmed present on the real rad
# archive.
_ACC_SUFFIX_RE = re.compile(r"_acc(\d+)h$", re.IGNORECASE)


def infer_deaccum_seconds(src_var, units):
    """
    Detect an ECMWF-style accumulated flux and return the number of
    seconds to divide by to convert it to an average irradiance/flux
    in W m-2, or None if the variable does not look accumulated.

    Two independent signals, checked together rather than alone:

      1. The variable name carries an explicit accumulation tag, e.g.
         "_acc1h" (seen on the real rad archive: "SSRD_GDS0_SFC_acc1h").
      2. The units string is dimensionally energy, not power:
         "J m**-2" (standard), or the nonstandard-but-real "W m**-2 s"
         seen here - dimensionally W*s = J, so "W m**-2 s" IS J m**-2
         despite the confusing "W" in the string. This is a genuine
         trap: the name says SSRD (radiation) and the units string
         starts with "W" (power), but the field is actually an
         accumulated energy over the period named in the variable.

    Requiring both signals (rather than units alone) avoids false
    positives on a file that is already instantaneous W m-2 but has a
    stray trailing character in its units string.
    """
    m = _ACC_SUFFIX_RE.search(src_var or "")
    hours = int(m.group(1)) if m else None

    u = (units or "").strip().lower()
    is_energy_units = u in ("j m**-2", "j m-2", "j/m^2", "j m^-2")
    is_mislabeled_power = u.startswith("w m") and u.endswith("s")

    if hours and (is_energy_units or is_mislabeled_power):
        return hours * 3600
    return None


class VarSpec:
    """How to find, read and label one output variable."""

    def __init__(self, name, subdir, path, src_var=None, prefer=None,
                 units=None, long_name=None, note=None):
        self.name = name
        self.subdir = subdir
        self.path = path
        self.src_var = src_var
        self.prefer = prefer or []
        self.units = units
        self.long_name = long_name
        self.note = note


def _p_2t(mm, dd, hh):
    return (f"2T/{mm:02d}/{dd:02d}/"
            f"2T.{MONTH_ABBR[mm-1]}{dd:02d}{hh:02d}.91-20clim.nc")


def _p_rad(mm, dd, hh):
    # NOTE: flat - no day subdirectory, unlike the other two.
    return f"rad/{mm:02d}/{mm:02d}{dd:02d}{hh:02d}.91-20radclim.nc"


def _p_ws10(mm, dd, hh):
    return f"windclimo/{mm:02d}/{dd:02d}/{mm:02d}{dd:02d}{hh:02d}.WS10.nc"


# The four families below arrive as tarballs whose members sit directly in
# a month directory ("01/STL1.jan1414.91-20clim.nc"), so extracting each
# family's tarballs with `tar -C <root>/<family>` yields exactly these
# paths. Names carry <DD><HH> (e.g. jan1414 = day 14, hour 14), the same
# convention as 2T; a wrong guess would show up as missing hours.

def _p_stl1(mm, dd, hh):
    return (f"STL1/{mm:02d}/"
            f"STL1.{MONTH_ABBR[mm-1]}{dd:02d}{hh:02d}.91-20clim.nc")


def _p_sd(mm, dd, hh):
    return f"sd/{mm:02d}/SD.{MONTH_ABBR[mm-1]}{dd:02d}{hh:02d}.91-20clim.nc"


def _p_skt(mm, dd, hh):
    return f"skt/{mm:02d}/SKT.{MONTH_ABBR[mm-1]}{dd:02d}{hh:02d}.91-20clim.nc"


def _p_wct(mm, dd, hh):
    # Same shape as rad: <MM><DD><HH> in the name, flat in the month dir.
    return f"wct/{mm:02d}/{mm:02d}{dd:02d}{hh:02d}.91-20wctclim.nc"


# Variables a plain build includes when --vars is not given. The newer
# families are opt-in (--vars STL1 ...), so an old command line keeps
# meaning what it always meant.
CORE_VARS = ("2T", "SSRD", "WS10")


def build_specs(rad_var=None):
    return {
        "2T": VarSpec(
            "2T", "2T", _p_2t, src_var="2T_GDS0_SFC",
            units="K", long_name="2 metre temperature"),
        "SSRD": VarSpec(
            "SSRD", "rad", _p_rad, src_var=rad_var,
            prefer=RAD_VAR_PREFERENCE,
            long_name="surface solar irradiance",
            note=("Source variable is accumulated energy (units "
                  "'W m**-2 s' = J m-2 per accumulation period); divide by "
                  "the period in seconds for average W m-2. Land-masked in "
                  "the source archive.")),
        "WS10": VarSpec(
            "WS10", "windclimo", _p_ws10, src_var=None,
            units="m s-1", long_name="10 metre wind speed",
            note="WS10 is precomputed speed; no u/v combination needed"),
        "STL1": VarSpec(
            "STL1", "STL1", _p_stl1, src_var=None,
            note=("Source variable is named STL1_GDS0_DBLY. Months not "
                  "delivered in the archive are NaN; see months_with_data.")),
        "SD": VarSpec(
            "SD", "sd", _p_sd, src_var=None,
            note=("Months not delivered in the archive are NaN; see "
                  "months_with_data.")),
        "SKT": VarSpec(
            "SKT", "skt", _p_skt, src_var=None,
            note=("Months not delivered in the archive are NaN; see "
                  "months_with_data.")),
        "WCT": VarSpec(
            "WCT", "wct", _p_wct, src_var=None,
            note=("The source files carry NO units and NO long_name, so "
                  "neither is known from the data. Verify the units "
                  "(check_wct_units.py) before relying on them.")),
    }


# ---------------------------------------------------------------------
# Discovery / IO
#
# Local root: plain os.listdir + xr.open_dataset(path). Fast, native,
# no extra dependency.
# HTTP root: directory index is HTML, parsed for href="NN/" entries;
# file reads go through fsspec so bytes are fetched as a file-like
# object rather than handing the raw URL to xr.open_dataset - the
# netCDF-C library treats a bare "http://" path as an OPeNDAP request
# and fails against a plain static file server (confirmed earlier in
# this project). That workaround is irrelevant for local paths, so the
# local branch skips it entirely.
# ---------------------------------------------------------------------

def available_months(root, subdir, fs=None):
    """Months published under <root>/<subdir>."""
    if is_url(root):
        with fs.open(f"{root}/{subdir}/", mode="rb") as f:
            html = f.read().decode("utf-8", errors="replace")
        return sorted({
            int(m) for m in re.findall(r'href="(?:.*/)?(\d{2})/"', html)
            if 1 <= int(m) <= 12
        })

    dirpath = os.path.join(root, subdir)
    if not os.path.isdir(dirpath):
        raise RuntimeError(f"Not a directory: {dirpath}")
    months = []
    for entry in os.listdir(dirpath):
        if (len(entry) == 2 and entry.isdigit()
                and os.path.isdir(os.path.join(dirpath, entry))):
            m = int(entry)
            if 1 <= m <= 12:
                months.append(m)
    return sorted(months)


def make_fs(root):
    """An fsspec HTTP filesystem for a URL root, or None for a local root."""
    if is_url(root):
        import fsspec
        return fsspec.filesystem("http")
    return None


# NetCDF-4 / HDF5 files must not be opened and read from several threads
# at once through the netCDF4 engine. xarray serialises the *reads* but
# gives every file its own lock for the *open*, and netCDF-C / HDF5 keep
# process-global state that is not safe under concurrent opens. The
# result is not an error but heap corruption: a segfault or
# "free(): invalid size", partway through a long run. Reproduced with
# files written as NetCDF-4 (a real WCT delivery was the trigger); classic
# (NetCDF-3) files, which never touch HDF5, survived the same load.
# So HDF5-format files take this lock for the whole open+read, while
# classic files stay fully parallel exactly as before.
_HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"
_HDF5_LOCK = threading.RLock()


def _is_hdf5_path(path):
    try:
        with open(path, "rb") as f:
            return f.read(8) == _HDF5_MAGIC
    except OSError:
        return False


@contextlib.contextmanager
def open_source(root, fs, spec, ts):
    """
    Open one hourly source file as an xarray Dataset, regardless of
    whether root is a local path or an http(s) URL. HDF5-format files are
    opened under a process-wide lock (see _HDF5_LOCK).
    """
    relpath = spec.path(ts.month, ts.day, ts.hour)
    if fs is not None:
        f = fs.open(f"{root}/{relpath}", mode="rb")
        try:
            hdf5 = f.read(8) == _HDF5_MAGIC
            f.seek(0)
            with (_HDF5_LOCK if hdf5 else contextlib.nullcontext()):
                with xr.open_dataset(f) as ds:
                    yield ds
        finally:
            f.close()
    else:
        path = os.path.join(root, relpath)
        with (_HDF5_LOCK if _is_hdf5_path(path) else contextlib.nullcontext()):
            with xr.open_dataset(path) as ds:
                yield ds


def _days_in_month(month, year=PLACEHOLDER_YEAR):
    nxt = dt.date(year + 1, 1, 1) if month == 12 else dt.date(year, month + 1, 1)
    return (nxt - dt.date(year, month, 1)).days


def build_times(months):
    return pd.DatetimeIndex([
        dt.datetime(PLACEHOLDER_YEAR, mm, dd, hh)
        for mm in months
        for dd in range(1, _days_in_month(mm) + 1)
        for hh in range(24)
    ])


def find_coord_names(ds):
    """
    Locate lat/lon coords without assuming g0_lat_0 / g0_lon_1 - the
    three archives were produced at different times and may not agree.
    """
    lat = lon = None
    for nm, v in ds.variables.items():
        u = str(v.attrs.get("units", "")).lower()
        if u.startswith("degrees_n"):
            lat = nm
        elif u.startswith("degrees_e"):
            lon = nm
    if lat is None or lon is None:
        for nm in ds.dims:
            low = str(nm).lower()
            if lat is None and "lat" in low:
                lat = nm
            if lon is None and "lon" in low:
                lon = nm
    if lat is None or lon is None:
        raise RuntimeError(
            f"Could not identify lat/lon coords. Variables: {list(ds.variables)}")
    return lat, lon


def pick_variable(ds, spec, latn, lonn):
    """
    Explicit spec.src_var wins; then spec.prefer in order; then, if
    exactly one 2-D (lat,lon) data variable exists, take it. If several
    exist and none matched, refuse rather than guess.
    """
    twod = [nm for nm, v in ds.data_vars.items()
            if latn in v.dims and lonn in v.dims]

    if spec.src_var:
        if spec.src_var not in ds.data_vars:
            raise RuntimeError(
                f"{spec.name}: variable '{spec.src_var}' not in file. "
                f"Available: {twod}")
        return spec.src_var

    for cand in spec.prefer:
        if cand in ds.data_vars:
            return cand

    if len(twod) == 1:
        return twod[0]

    raise RuntimeError(
        f"{spec.name}: cannot choose a variable automatically. Candidates: "
        f"{twod}. Re-run with an explicit name (e.g. --rad-var <NAME>).")


def inspect(root, spec, fs=None):
    """Open the first readable file for this variable and report contents."""
    months = available_months(root, spec.subdir, fs)
    for ts in build_times(months):
        try:
            with open_source(root, fs, spec, ts) as ds:
                latn, lonn = find_coord_names(ds)
                return {
                    "months": months,
                    "lat_name": latn, "lon_name": lonn,
                    "nlat": ds.sizes[latn], "nlon": ds.sizes[lonn],
                    "lats": ds[latn].values.astype("float32"),
                    "lons": ds[lonn].values.astype("float32"),
                    "data_vars": {
                        nm: dict(units=v.attrs.get("units"),
                                 long_name=v.attrs.get("long_name"))
                        for nm, v in ds.data_vars.items()},
                    "sample_ts": ts,
                    "sample_path": f"{root}/{spec.path(ts.month, ts.day, ts.hour)}",
                }
        except Exception:
            continue
    raise RuntimeError(f"{spec.name}: no readable source file found under {root}.")


def list_vars(names=None, rad_var=None, root=DEFAULT_ROOT):
    """Print what is actually inside each archive's files."""
    fs = make_fs(root)
    specs = build_specs(rad_var)
    for nm in (names or specs):
        spec = specs[nm]
        print(f"\n=== {nm}  ({spec.subdir}) ===")
        try:
            info = inspect(root, spec, fs)
        except Exception as exc:
            print(f"  FAILED: {exc}")
            continue
        print(f"  sample: {info['sample_path']}")
        print(f"  months: {info['months']}")
        print(f"  grid:   {info['nlat']} x {info['nlon']} "
              f"({info['lat_name']}, {info['lon_name']})")
        print(f"  lat:    {info['lats'][0]} .. {info['lats'][-1]}")
        print(f"  lon:    {info['lons'][0]} .. {info['lons'][-1]}")
        print("  data variables:")
        for v, meta in info["data_vars"].items():
            print(f"    - {v:26s} {str(meta['units'] or ''):12s} "
                  f"{meta['long_name'] or ''}")


# ---------------------------------------------------------------------
# Conversion
# ---------------------------------------------------------------------

def _fetch(root, fs, spec, src_var, ts):
    """One hourly global field, or None if absent/unreadable."""
    try:
        with open_source(root, fs, spec, ts) as ds:
            return np.asarray(ds[src_var].values, dtype="float32")
    except Exception:
        return None


def _resolve_source_var(root, fs, spec, sample_ts):
    """(variable name, units, long_name) read from one real source file."""
    with open_source(root, fs, spec, sample_ts) as ds:
        latn, lonn = find_coord_names(ds)
        sv = pick_variable(ds, spec, latn, lonn)
        return sv, ds[sv].attrs.get("units"), ds[sv].attrs.get("long_name")


def _grid_group_name(lats):
    """Self-describing group label from grid spacing, e.g. 'grid_0p25'."""
    return "grid_" + f"{abs(float(lats[1] - lats[0])):.2f}".replace(".", "p")


def _time_values(times):
    return ((times - pd.Timestamp(f"{PLACEHOLDER_YEAR}-01-01"))
            // pd.Timedelta("1h")).values.astype("int64")


def _write_group_coords(grp, lats, lons, tvals, chunks):
    """lat / lon / time coordinate arrays for one grid group."""
    nla, nlo, nt = len(lats), len(lons), len(tvals)
    ct, cla, clo = chunks

    # fill_value=None on coords is essential: zarr v2 defaults it to 0
    # and xarray then reads 0 as missing - which would silently NaN out
    # lon=0.0 (the prime meridian IS on this grid) and time index 0.
    zlat = _make_array(grp, "lat", (nla,), (nla,), "float32", fill_value=None)
    zlat[:] = lats
    zlat.attrs["_ARRAY_DIMENSIONS"] = ["lat"]
    zlat.attrs["units"] = "degrees_north"

    zlon = _make_array(grp, "lon", (nlo,), (nlo,), "float32", fill_value=None)
    zlon[:] = lons
    zlon.attrs["_ARRAY_DIMENSIONS"] = ["lon"]
    zlon.attrs["units"] = "degrees_east"

    ztime = _make_array(grp, "time", (nt,), (nt,), "int64", fill_value=None)
    ztime[:] = tvals
    ztime.attrs["_ARRAY_DIMENSIONS"] = ["time"]
    ztime.attrs["units"] = f"hours since {PLACEHOLDER_YEAR}-01-01 00:00:00"
    ztime.attrs["calendar"] = "standard"
    grp.attrs["chunking_rationale"] = (
        f"chunked ({ct},{cla},{clo}) so a point timeseries costs "
        f"{int(np.ceil(nt/ct))} request(s) instead of {nt}")


def _apply_meta(arr, attrs):
    """
    Set operator-supplied attributes. Supplying units records who said so
    (units_source) and clears the "units unknown" note - unless the
    operator supplied units_source themselves, which then wins regardless
    of argument order (e.g. "inferred from data, not confirmed").
    """
    attrs = dict(attrs or {})
    for k, v in attrs.items():
        arr.attrs[k] = v
    if "units" in attrs:
        if "units_source" not in attrs:
            arr.attrs["units_source"] = "supplied by the operator (--meta)"
        if "units_note" in arr.attrs:
            del arr.attrs["units_note"]


def _define_array(grp, nm, spec, root, sv, src_units, src_long, shape, chunks,
                  compressor, deaccumulate, meta, verbose):
    """
    Create one data array with its full attribute set. Shared by build
    and append so the two can never disagree about what a variable's
    metadata looks like. Returns (array, deaccumulation_seconds_or_None).
    """
    deaccum_seconds = (infer_deaccum_seconds(sv, src_units)
                       if deaccumulate else None)
    out_units = spec.units or src_units or ""
    if deaccum_seconds:
        out_units = "W m-2"
        if verbose:
            print(f"  {nm} <- '{sv}' (units: {src_units}) "
                  f"** ACCUMULATED FLUX DETECTED: dividing by "
                  f"{deaccum_seconds}s -> stored as W m-2 **")
    elif verbose:
        print(f"  {nm} <- '{sv}' (units: {src_units})")

    nt, nla, nlo = shape
    a = _make_array(grp, nm, shape, chunks, "float32",
                    compressor=compressor, fill_value=np.nan)
    a.attrs["_ARRAY_DIMENSIONS"] = ["time", "lat", "lon"]
    a.attrs["units"] = out_units
    a.attrs["long_name"] = spec.long_name or src_long or nm
    a.attrs["source_variable"] = sv
    a.attrs["source_units"] = src_units or ""
    # A local root is recorded as an absolute path: a relative one is
    # meaningless in a store's metadata once you are in another directory.
    a.attrs["source_layout"] = (
        f"{root if is_url(root) else os.path.abspath(root)}/"
        f"{spec.path(1, 1, 0)}")
    a.attrs["grid"] = f"{nla}x{nlo}"
    # Record the *option*, not just the outcome, so a verifier can tell
    # "deliberately left raw" from "should have been converted but wasn't".
    a.attrs["deaccumulate_option"] = bool(deaccumulate)
    if deaccum_seconds:
        a.attrs["deaccumulated"] = True
        a.attrs["deaccumulation_seconds"] = deaccum_seconds
        a.attrs["note"] = (
            (spec.note + " " if spec.note else "") +
            f"Converted from accumulated '{src_units}' to average "
            f"W m-2 by dividing by {deaccum_seconds}s. Set "
            f"deaccumulate=False to store raw accumulated values.")
    elif spec.note:
        a.attrs["note"] = spec.note

    # A variable whose files declare no units must not look like it has
    # them: record the gap rather than leaving a bare empty string.
    if not out_units:
        a.attrs["units_note"] = ("source files carry no units attribute; "
                                 "units are unknown")
    _apply_meta(a, meta)
    return a, deaccum_seconds


def _fill_array(arr, root, fs, spec, sv, times, deaccum_seconds,
                quantize_digits, max_workers, verbose=True, label=""):
    """
    Read every hour for one variable and write it, one time-chunk at a
    time (so peak memory is one chunk of the full grid). Missing hours
    stay NaN. Returns {"missing": n, "months_with_data": [...]}.
    """
    nt = len(times)
    nla, nlo = arr.shape[1], arr.shape[2]
    step = arr.chunks[0]
    got = collections.Counter()
    missing = 0
    for start in range(0, nt, step):
        stop = min(start + step, nt)
        block = times[start:stop]
        buf = np.full((stop - start, nla, nlo), np.nan, dtype="float32")
        with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
            futs = {ex.submit(_fetch, root, fs, spec, sv, ts): i
                    for i, ts in enumerate(block)}
            for fut in cf.as_completed(futs):
                i = futs[fut]
                field = fut.result()
                if field is None:
                    missing += 1
                else:
                    buf[i] = field
                    got[int(block[i].month)] += 1
        if deaccum_seconds:
            buf = buf / deaccum_seconds
        if quantize_digits is not None:
            buf = np.round(buf, quantize_digits)
        arr[start:stop, :, :] = buf
        if verbose:
            print(f"  {label or spec.name}: wrote {stop}/{nt} timesteps "
                  f"({100*stop/nt:.0f}%)")
    return {"missing": missing, "months_with_data": sorted(got)}


def _record_coverage(arr, nt, stats):
    """Store which months actually have data, so gaps are visible later."""
    arr.attrs["months_with_data"] = [int(m) for m in stats["months_with_data"]]
    arr.attrs["missing_hours"] = int(stats["missing"])
    arr.attrs["expected_hours"] = int(nt)


def convert(out_path, var_names=None, rad_var=None, root=DEFAULT_ROOT,
            chunks=DEFAULT_CHUNKS, max_workers=16, clevel=5,
            quantize_digits=None, zarr_format=2, deaccumulate=True,
            meta=None, verbose=True):
    """
    Build one Zarr store containing every requested variable.

    root : local filesystem path (default: the Qumulo mount at
        /qumulo/azcot) or an http(s) URL - auto-detected.

    var_names : defaults to CORE_VARS (2T, SSRD, WS10). The newer
        families (STL1, SD, SKT, WCT) are opt-in. To add them to a store
        that already exists, use append() instead - this function
        creates the output from scratch and overwrites what is there.

    All variables on one grid share a group; differing grids get their
    own group and are never regridded (resampling climatology fields is
    a scientific decision, not a plumbing one).

    deaccumulate : if True (default), a variable whose name carries an
        ECMWF-style accumulation tag ("..._acc1h") AND whose units are
        dimensionally energy (J m-2, or the confusing-but-real
        "W m**-2 s" = J m-2) is divided by the accumulation period in
        seconds, converting it to an average flux in W m-2. Confirmed
        needed for the real rad/SSRD archive: its variable is literally
        named "SSRD_GDS0_SFC_acc1h" with units "W m**-2 s" - i.e. it is
        accumulated energy, not instantaneous irradiance, despite the
        "W" in the units string. Left unconverted, values are 3600x
        too large for a 1h accumulation and would not look obviously
        wrong. Set False to store the raw accumulated values as-is.

    meta : {variable: {attr: value}} overrides applied to array
        attributes, e.g. {"WCT": {"units": "degF"}}.
    """
    specs = build_specs(rad_var)
    var_names = list(var_names) if var_names else list(CORE_VARS)
    meta = meta or {}
    fs = make_fs(root)

    if verbose:
        print(f"Reading from: {root} ({'HTTP' if fs else 'local filesystem'})")
        print("Inspecting sources...")
    info = {nm: inspect(root, specs[nm], fs) for nm in var_names}
    if verbose:
        for nm in var_names:
            print(f"  {nm}: grid {info[nm]['nlat']}x{info[nm]['nlon']}, "
                  f"months {info[nm]['months']}")

    # Group variables by grid. The archives genuinely differ: 2T and
    # WS10 are 0.25 deg (121x1440) while rad/SSRD is 0.10 deg
    # (301x3600) over the same lat 90..60, lon -180..180 domain.
    #
    # They are NOT regridded onto a common grid, deliberately:
    #   * 0.25/0.10 = 2.5 is not an integer ratio, so points coincide
    #     only every 0.5 deg - any merge is real interpolation.
    #   * Downsampling SSRD would destroy resolution the provider paid
    #     to produce; upsampling 2T/WS10 would inflate them 6.25x for
    #     zero new information.
    #   * Regridding a flux like solar irradiance properly wants
    #     area-weighted/conservative remapping, which is a scientific
    #     choice for the data owner, not something to do implicitly.
    #   * For point timeseries - the whole purpose here - a common grid
    #     is unnecessary: pick the nearest cell in each grid.
    #
    # So each distinct grid becomes its own group in one store.
    grids = {}
    for nm in var_names:
        key = (info[nm]["nlat"], info[nm]["nlon"])
        grids.setdefault(key, []).append(nm)

    info_lats = {(info[nm]["nlat"], info[nm]["nlon"]): info[nm]["lats"]
                 for nm in var_names}

    def group_name(nlat, nlon):
        return _grid_group_name(info_lats[(nlat, nlon)])

    single_grid = len(grids) == 1

    if verbose:
        print()
        if single_grid:
            print("All variables share one grid -> writing flat at store root.")
        else:
            print(f"{len(grids)} distinct grids -> one group each "
                  f"(no regridding):")
            for (nla, nlo), members in grids.items():
                dlat = abs(float(info_lats[(nla, nlo)][1]
                                 - info_lats[(nla, nlo)][0]))
                print(f"  {group_name(nla, nlo):12s} {nla}x{nlo} "
                      f"(~{dlat:.2f} deg): {', '.join(members)}")

    # Time axis: union of published months across all requested
    # variables, so a variable missing a month becomes NaN rather than
    # shifting the time axis for everything else. Written identically
    # into every group so the groups stay time-aligned.
    months = sorted(set().union(*(set(info[nm]["months"]) for nm in var_names)))
    times = build_times(months)
    nt = len(times)
    tvals = _time_values(times)

    ct_req, cla, clo = chunks
    ct = min(ct_req, nt)

    if verbose:
        print(f"\nTime axis: {nt} hourly steps over months {months}")
        for (nla, nlo), members in grids.items():
            nchunk = (int(np.ceil(nt/ct)) * int(np.ceil(nla/cla))
                      * int(np.ceil(nlo/clo)))
            raw = nt*nla*nlo*4/1e9
            print(f"  {nla}x{nlo}: {raw:.2f} GB raw/variable "
                  f"(~{raw/1.8:.2f} GB zstd), {nchunk:,} objects/variable")
        print(f"Chunks ({ct},{cla},{clo}) ~ {ct*cla*clo*4/1e6:.2f} MB")
        print(f"Point-year cost: {int(np.ceil(nt/ct))} request(s) vs {nt} today")
        print(f"Peak RAM per block ~ up to "
              f"{ct*max(k[0]*k[1] for k in grids)*4/1e9:.2f} GB")

    compressor = _make_compressor(clevel, zarr_format)
    zroot = _open_group(out_path, zarr_format)
    zroot.attrs["title"] = (
        "azcot 1991-2020 hourly climatology, timeseries-optimized")
    zroot.attrs["comment"] = (
        "The year in the time coordinate is a placeholder and carries no "
        "meaning; only month/day/hour are real. NaN = hour not published.")
    if not single_grid:
        zroot.attrs["layout"] = (
            "One group per source grid; grids are NOT interpolated onto a "
            "common mesh. Open a group with "
            "xr.open_zarr(url, group='grid_0p25').")
    zroot.attrs["groups"] = {group_name(*k): {"shape": [nt, k[0], k[1]],
                                              "variables": v}
                             for k, v in grids.items()}

    made = {}
    for (nla, nlo), members in grids.items():
        grp = zroot if single_grid else zroot.create_group(group_name(nla, nlo))
        _write_group_coords(grp, info_lats[(nla, nlo)],
                            info[members[0]]["lons"], tvals, (ct, cla, clo))
        for nm in members:
            spec = specs[nm]
            sv, src_units, src_long = _resolve_source_var(
                root, fs, spec, info[nm]["sample_ts"])
            a, dsec = _define_array(
                grp, nm, spec, root, sv, src_units, src_long,
                (nt, nla, nlo), (ct, cla, clo), compressor, deaccumulate,
                meta.get(nm), verbose)
            made[nm] = (a, sv, dsec)

    stats_by_var = {}
    for nm in var_names:
        a, sv, dsec = made[nm]
        stats = _fill_array(a, root, fs, specs[nm], sv, times, dsec,
                            quantize_digits, max_workers, verbose, nm)
        _record_coverage(a, nt, stats)
        stats_by_var[nm] = stats

    zarr.consolidate_metadata(out_path)   # 1 metadata GET, not thousands

    if verbose:
        print(f"\nDone -> {out_path}")
        for nm in var_names:
            dsec = made[nm][2]
            tag = f" (deaccumulated /{dsec}s -> W m-2)" if dsec else ""
            st = stats_by_var[nm]
            print(f"  {nm}: {st['missing']} missing/unreadable hours, "
                  f"months with data {st['months_with_data']}{tag}")
        if not single_grid:
            print("\nRead with:")
            for (nla, nlo), members in grids.items():
                print(f"  xr.open_zarr(url, group='{group_name(nla, nlo)}')"
                      f"  # {', '.join(members)}")
    return out_path


# ---------------------------------------------------------------------
# Adding to / annotating an existing store
# ---------------------------------------------------------------------

def _open_for_append(path):
    """
    Open an existing store read-write. zarr 3 would otherwise trust the
    consolidated metadata, which goes stale the moment an array is
    added; it is regenerated at the end of every append.
    """
    if ZARR2:
        return zarr.open_group(path, mode="r+")
    return zarr.open_group(path, mode="r+", use_consolidated=False)


def _read_store_times(zroot, groups):
    """The store's own time axis is the authority; refuse one we didn't make."""
    first = next(iter(groups))
    tv = np.asarray(zroot[first]["time"][:], dtype="int64")
    times = pd.DatetimeIndex(pd.Timestamp(f"{PLACEHOLDER_YEAR}-01-01")
                             + pd.to_timedelta(tv, unit="h"))
    store_months = sorted(set(int(m) for m in times.month))
    if not times.equals(build_times(store_months)):
        raise RuntimeError(
            "the store's time axis is not the standard hourly axis for its "
            "months; refusing to modify a store whose axis this tool did "
            "not make")
    return times, store_months


def _store_zarr_format(zroot):
    return getattr(getattr(zroot, "metadata", None), "zarr_format", 2)


def _adoptable(arr, sv, shape, deaccumulate):
    """
    Is this unregistered array a COMPLETE result of an earlier interrupted
    run that matches what is being asked for now? months_with_data is
    written only after a variable's last block has been stored, so its
    presence is the completion marker; a half-written array lacks it.
    """
    a = arr.attrs
    return ("months_with_data" in a
            and a.get("source_variable") == sv
            and tuple(arr.shape) == tuple(shape)
            and bool(a.get("deaccumulate_option")) == bool(deaccumulate))


def append(out_path, var_names, root=DEFAULT_ROOT, rad_var=None, replace=(),
           meta=None, chunks=DEFAULT_CHUNKS, max_workers=16, clevel=5,
           quantize_digits=None, deaccumulate=True, resume=False,
           verbose=True):
    """
    Add variables to an EXISTING store without touching what is in it.

    Nothing already in the store is rewritten. Each new variable is
    placed in the existing group whose lat/lon it matches exactly, or in
    a new group if none does (never regridded). It must fit the store's
    time axis: months in the source that the store does not have are
    refused rather than silently dropped, while months the store has but
    the source lacks are left NaN and recorded in months_with_data.

    replace : names of variables that may be overwritten. Anything else
        already in the store is refused, so a repeated run cannot clobber
        data.

    The root attribute 'groups' and the consolidated metadata are only
    updated once every new variable has been written, so an interrupted
    run leaves unregistered (invisible) arrays rather than a store that
    advertises variables it does not fully contain. (Variables being
    replaced are un-registered first, before their data is deleted.)

    Re-running after such an interruption: unregistered leftovers of the
    variables being added are deleted and rewritten (they were never part
    of the store). With resume=True, a leftover that was written
    completely and matches the request (same source variable, shape and
    deaccumulation option) is kept instead of re-read; a half-written one
    is always deleted.
    """
    if not var_names:
        raise ValueError("append() needs an explicit list of variables")
    specs = build_specs(rad_var)
    unknown = [v for v in var_names if v not in specs]
    if unknown:
        raise ValueError(f"unknown variable(s) {unknown}; known: {sorted(specs)}")
    meta = meta or {}
    fs = make_fs(root)

    zroot = _open_for_append(out_path)
    raw_groups = dict(zroot.attrs.get("groups") or {})
    if not raw_groups:
        raise RuntimeError(
            f"{out_path} has no 'groups' attribute, so it is not a grouped "
            f"store this tool can append to.")
    groups = {g: {"shape": list(m["shape"]), "variables": list(m["variables"])}
              for g, m in raw_groups.items()}
    existing = {v: g for g, m in groups.items() for v in m["variables"]}

    clash = [v for v in var_names if v in existing and v not in replace]
    if clash:
        raise RuntimeError(
            f"already in the store: {clash}. Nothing was changed. Pass "
            f"--replace {' '.join(clash)} to overwrite deliberately.")

    times, store_months = _read_store_times(zroot, groups)
    nt = len(times)
    tvals = _time_values(times)

    if verbose:
        print(f"Appending to {out_path}")
        print(f"  existing groups: "
              f"{ {g: m['variables'] for g, m in groups.items()} }")
        print(f"  time axis: {nt} steps, months {store_months}")
        print(f"Reading from: {root} ({'HTTP' if fs else 'local filesystem'})")

    info = {nm: inspect(root, specs[nm], fs) for nm in var_names}
    for nm in var_names:
        extra = sorted(set(info[nm]["months"]) - set(store_months))
        if extra:
            raise RuntimeError(
                f"{nm}: source has month(s) {extra} that are not on the "
                f"store's time axis {store_months}. Nothing was changed. "
                f"Adding them would need a rebuild with a wider axis.")
        lack = sorted(set(store_months) - set(info[nm]["months"]))
        if verbose:
            print(f"  {nm}: grid {info[nm]['nlat']}x{info[nm]['nlon']}, "
                  f"source months {info[nm]['months']}"
                  + (f"  -> months {lack} will be NaN" if lack else ""))

    def match_group(nm):
        la, lo = info[nm]["lats"], info[nm]["lons"]
        for g in groups:
            glat = np.asarray(zroot[g]["lat"][:], dtype="float64")
            glon = np.asarray(zroot[g]["lon"][:], dtype="float64")
            if (len(glat) == len(la) and len(glon) == len(lo)
                    and np.allclose(glat, la, atol=1e-3)
                    and np.allclose(glon, lo, atol=1e-3)):
                return g
        return None

    plan, pending = {}, []      # pending: new groups to create
    for nm in var_names:
        g = match_group(nm)
        if g is None:
            la, lo = info[nm]["lats"], info[nm]["lons"]
            for pg in pending:
                if (len(pg["lats"]) == len(la) and len(pg["lons"]) == len(lo)
                        and np.allclose(pg["lats"], la, atol=1e-3)
                        and np.allclose(pg["lons"], lo, atol=1e-3)):
                    g = pg["name"]
                    break
        if g is None:
            la, lo = info[nm]["lats"], info[nm]["lons"]
            base = _grid_group_name(la)
            taken = set(groups) | {pg["name"] for pg in pending}
            # Same spacing but different shape/domain must not collide
            # with (or silently reuse) an existing group of that name.
            name = base if base not in taken else f"{base}_{len(la)}x{len(lo)}"
            pending.append({"name": name, "lats": la, "lons": lo})
            g = name
        plan[nm] = g

    compressor = _make_compressor(clevel, _store_zarr_format(zroot))

    # Deliberate overwrites. ORDER MATTERS: stop advertising the variable
    # first, then delete its data, then refresh the consolidated metadata.
    # Deleting first (as an earlier version did) meant an interruption
    # during the long rewrite left the store listing a variable whose data
    # was half destroyed. This way an interruption leaves the variable
    # simply absent, which a re-run restores. The variable is therefore
    # unavailable to readers while it is being rewritten.
    replaced = [nm for nm in var_names if nm in existing]
    if replaced:
        for nm in replaced:
            groups[existing[nm]]["variables"].remove(nm)
        zroot.attrs["groups"] = groups            # 1. no longer advertised
        for nm in replaced:
            del zroot[existing[nm]][nm]           # 2. now remove the data
        zarr.consolidate_metadata(out_path)       # 3. readers stop seeing it
        if verbose:
            for nm in replaced:
                print(f"  {nm}: existing array removed (--replace); absent "
                      f"from the store until the rewrite completes")

    for pg in pending:
        if pg["name"] in zroot:            # unregistered group from a crashed run
            del zroot[pg["name"]]
            if verbose:
                print(f"  group {pg['name']}: removing unregistered leftover "
                      f"from an interrupted run")
        grp = zroot.create_group(pg["name"])
        ct = min(chunks[0], nt)
        _write_group_coords(grp, pg["lats"], pg["lons"], tvals,
                            (ct, chunks[1], chunks[2]))
        groups[pg["name"]] = {"shape": [nt, len(pg["lats"]), len(pg["lons"])],
                              "variables": []}
        if verbose:
            print(f"  new group {pg['name']} "
                  f"({len(pg['lats'])}x{len(pg['lons'])})")

    def chunks_for(g):
        for v in groups[g]["variables"]:
            try:
                return tuple(int(c) for c in zroot[g][v].chunks)
            except Exception:
                pass
        return (min(chunks[0], nt), chunks[1], chunks[2])

    results = {}
    for nm in var_names:
        g = plan[nm]
        grp = zroot[g]
        spec = specs[nm]
        nla, nlo = info[nm]["nlat"], info[nm]["nlon"]
        ch = chunks_for(g)
        sv, src_units, src_long = _resolve_source_var(
            root, fs, spec, info[nm]["sample_ts"])

        # An array on disk that the registry does not list was left by an
        # interrupted run (a registered one was refused or replaced above).
        if nm in set(grp.array_keys()) and nm not in groups[g]["variables"]:
            if resume and _adoptable(grp[nm], sv, (nt, nla, nlo), deaccumulate):
                arr = grp[nm]
                _apply_meta(arr, meta.get(nm))
                groups[g]["variables"].append(nm)
                dsec = (arr.attrs.get("deaccumulation_seconds")
                        if arr.attrs.get("deaccumulated") else None)
                results[nm] = (g, {
                    "missing": int(arr.attrs["missing_hours"]),
                    "months_with_data": list(arr.attrs["months_with_data"])},
                    dsec)
                if verbose:
                    print(f"  {nm}: keeping the complete array an interrupted "
                          f"run wrote (--resume); not re-reading it")
                continue
            del grp[nm]
            if verbose:
                print(f"  {nm}: removing unregistered leftover from an "
                      f"interrupted run; rewriting it")

        a, dsec = _define_array(
            grp, nm, spec, root, sv, src_units, src_long, (nt, nla, nlo), ch,
            compressor, deaccumulate, meta.get(nm), verbose)
        if verbose:
            print(f"  {nm} -> group {g}, chunks {ch}")
        stats = _fill_array(a, root, fs, spec, sv, times, dsec,
                            quantize_digits, max_workers, verbose, nm)
        _record_coverage(a, nt, stats)
        groups[g]["variables"].append(nm)
        results[nm] = (g, stats, dsec)

    zroot.attrs["groups"] = groups              # register only when complete
    zarr.consolidate_metadata(out_path)

    if verbose:
        print(f"\nDone. Appended to {out_path}")
        for nm, (g, st, dsec) in results.items():
            tag = f" (deaccumulated /{dsec}s)" if dsec else ""
            print(f"  {nm} in {g}: {st['missing']} missing hours, "
                  f"months with data {st['months_with_data']}{tag}")
    return results


def fill(out_path, var_names, root=DEFAULT_ROOT, rad_var=None, months=None,
         max_workers=16, quantize_digits=None, verbose=True):
    """
    Write months that were NaN because their source files had not arrived,
    now that they have - without touching anything else.

    For each variable, the months to fill are those whose source files now
    exist and that the array's months_with_data attribute does not list
    (or exactly `months`, if given). Only those hours are read and written,
    in place. Nothing is deleted, so the array is valid at every moment:
    an interruption leaves it with the months it had, and re-running
    fills the rest. months_with_data is updated only after a month is
    written, and a month not listed there has no data in the store, so
    overwriting its whole slice can never destroy anything.

    Each month is checked against what the array already holds before
    it is written: the source variable's units must equal the recorded
    source_units, and the grid must equal the group's. A later delivery
    that quietly changed either is refused, not stored.
    """
    if not var_names:
        raise ValueError("fill() needs an explicit list of variables")
    specs = build_specs(rad_var)
    fs = make_fs(root)
    zroot = _open_for_append(out_path)
    groups = dict(zroot.attrs.get("groups") or {})
    if not groups:
        raise RuntimeError(f"{out_path} has no 'groups' attribute")
    where = {v: g for g, m in groups.items() for v in m["variables"]}
    absent = [v for v in var_names if v not in where]
    if absent:
        raise RuntimeError(f"not in the store: {absent}. Use --append to add "
                           f"a variable; --fill only completes existing ones.")
    times, store_months = _read_store_times(zroot, groups)
    month_of = np.asarray(times.month)

    summary = {}
    for nm in var_names:
        spec = specs[nm]
        g = where[nm]
        arr = zroot[g][nm]
        have = arr.attrs.get("months_with_data")
        if have is None and months is None:
            raise RuntimeError(
                f"{nm} has no months_with_data attribute (it predates it), so "
                f"which months are empty is not recorded. Pass --months "
                f"explicitly for months you know are empty.")
        have = set(int(m) for m in (have or []))
        info = inspect(root, spec, fs)
        src_months = set(info["months"])

        if months is None:
            todo = sorted((src_months & set(store_months)) - have)
        else:
            todo = sorted(set(int(m) for m in months))
            bad = [m for m in todo if m not in store_months or m not in src_months]
            if bad:
                raise RuntimeError(
                    f"{nm}: month(s) {bad} are not both on the store's time "
                    f"axis {store_months} and present in the source "
                    f"{sorted(src_months)}")
            held = [m for m in todo if m in have]
            if held:
                raise RuntimeError(
                    f"{nm}: month(s) {held} already have data. --fill only "
                    f"completes empty months; to refresh data use --replace.")
        if not todo:
            if verbose:
                print(f"  {nm}: nothing to fill (months with data "
                      f"{sorted(have)}; source months {sorted(src_months)})")
            summary[nm] = []
            continue

        sv = arr.attrs["source_variable"]
        want_units = arr.attrs.get("source_units", "")
        dsec = (arr.attrs.get("deaccumulation_seconds")
                if arr.attrs.get("deaccumulated") else None)
        glat = np.asarray(zroot[g]["lat"][:], dtype="float64")
        glon = np.asarray(zroot[g]["lon"][:], dtype="float64")
        nla, nlo = arr.shape[1], arr.shape[2]

        written, total_read = [], 0
        for m in todo:
            idx = np.flatnonzero(month_of == m)
            i0, i1 = int(idx[0]), int(idx[-1]) + 1
            if i1 - i0 != len(idx):
                raise RuntimeError(f"month {m} is not contiguous on the time axis")
            block = times[i0:i1]

            # Guard: a month delivered later must match what is stored.
            sample = None
            for ts in block:
                try:
                    with open_source(root, fs, spec, ts) as ds:
                        latn, lonn = find_coord_names(ds)
                        got_units = ds[sv].attrs.get("units") or ""
                        la = np.asarray(ds[latn].values, dtype="float64")
                        lo = np.asarray(ds[lonn].values, dtype="float64")
                    sample = ts
                    break
                except Exception:
                    continue
            if sample is None:
                if verbose:
                    print(f"  {nm}: month {m:02d}: no readable file, skipped")
                continue
            if got_units != want_units:
                raise RuntimeError(
                    f"{nm}: month {m:02d} files declare units {got_units!r} but "
                    f"the stored data came from units {want_units!r}. Refusing "
                    f"to mix them. Nothing from this month was written.")
            if (la.shape != glat.shape or lo.shape != glon.shape
                    or not np.allclose(la, glat, atol=1e-3)
                    or not np.allclose(lo, glon, atol=1e-3)):
                raise RuntimeError(
                    f"{nm}: month {m:02d} files are on a different grid from "
                    f"the stored data. Nothing from this month was written.")

            buf = np.full((i1 - i0, nla, nlo), np.nan, dtype="float32")
            n_ok = 0
            with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
                futs = {ex.submit(_fetch, root, fs, spec, sv, ts): i
                        for i, ts in enumerate(block)}
                for fut in cf.as_completed(futs):
                    field = fut.result()
                    if field is not None:
                        buf[futs[fut]] = field
                        n_ok += 1
            if n_ok == 0:
                if verbose:
                    print(f"  {nm}: month {m:02d}: no readable hours, skipped")
                continue
            if dsec:
                buf = buf / dsec
            if quantize_digits is not None:
                buf = np.round(buf, quantize_digits)
            arr[i0:i1, :, :] = buf           # whole month was empty: safe
            written.append(m)
            total_read += n_ok
            if verbose:
                print(f"  {nm}: month {m:02d}: wrote {n_ok}/{i1 - i0} hours")

        if written:
            arr.attrs["months_with_data"] = sorted(have | set(written))
            old_missing = arr.attrs.get("missing_hours")
            if old_missing is not None:
                arr.attrs["missing_hours"] = max(0, int(old_missing) - total_read)
        summary[nm] = written
        # Refresh after every variable, not only at the end, so the
        # published metadata never lags far behind the data if a later
        # variable is interrupted.
        zarr.consolidate_metadata(out_path)

    if verbose:
        print(f"\nDone. Filled: {summary}")
    return summary


def annotate(out_path, meta, verbose=True):
    """
    Change attributes on arrays in an existing store, writing no data.
    Used e.g. to record units once they have been established.
    """
    if not meta:
        raise ValueError("nothing to annotate")
    zroot = _open_for_append(out_path)
    groups = dict(zroot.attrs.get("groups") or {})
    where = {v: g for g, m in groups.items() for v in m["variables"]}
    missing = [v for v in meta if v not in where]
    if missing:
        raise RuntimeError(f"not in the store: {missing}; have {sorted(where)}")
    for nm, attrs in meta.items():
        arr = zroot[where[nm]][nm]
        _apply_meta(arr, attrs)
        if verbose:
            for k, v in attrs.items():
                print(f"  {nm}.{k} = {v!r}")
    zarr.consolidate_metadata(out_path)
    if verbose:
        print("metadata re-consolidated")


def parse_meta(items, valid):
    """['WCT.units=degF', ...] -> {'WCT': {'units': 'degF'}}"""
    out = {}
    for it in items or []:
        key, sep, val = it.partition("=")
        var, dot, attr = key.partition(".")
        if not (sep and dot and var and attr):
            raise SystemExit(f"--meta expects VAR.attr=value, got {it!r}")
        if var not in valid:
            raise SystemExit(f"--meta: unknown variable {var!r}; "
                             f"known: {sorted(valid)}")
        out.setdefault(var, {})[attr] = val
    return out


def main():
    # A crash inside a C library (netCDF-C, HDF5, blosc) is not a Python
    # exception and normally leaves only "Segmentation fault (core
    # dumped)". faulthandler prints every thread's Python stack at that
    # moment, which shows which function and file it died in.
    import faulthandler
    faulthandler.enable()
    known = sorted(build_specs())
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", help="output .zarr path")
    p.add_argument("--root", default=DEFAULT_ROOT,
                   help=f"source root: local path or http(s) URL "
                        f"(default: {DEFAULT_ROOT})")
    p.add_argument("--vars", nargs="*", default=None,
                   choices=known + ["all"],
                   help=f"variables to include. Build default: "
                        f"{' '.join(CORE_VARS)}. 'all' = every known "
                        f"variable. Required with --append.")
    p.add_argument("--append", action="store_true",
                   help="add --vars to the EXISTING store at --out instead "
                        "of building one (existing data is not touched)")
    p.add_argument("--replace", nargs="*", default=[], choices=known,
                   help="with --append: variables allowed to be overwritten")
    p.add_argument("--fill", action="store_true",
                   help="complete existing variables with months whose source "
                        "files have arrived since (writes only those months; "
                        "deletes nothing)")
    p.add_argument("--months", type=int, nargs="*", default=None,
                   help="with --fill: fill exactly these months instead of "
                        "detecting them")
    p.add_argument("--resume", action="store_true",
                   help="with --append: after an interrupted run, keep "
                        "variables it wrote completely instead of re-reading "
                        "them (half-written ones are always redone)")
    p.add_argument("--annotate", action="store_true",
                   help="only change array attributes in the existing store "
                        "(via --meta); writes no data")
    p.add_argument("--meta", action="append", default=[],
                   metavar="VAR.attr=value",
                   help="set an array attribute, e.g. WCT.units=degF; "
                        "repeatable")
    p.add_argument("--rad-var", default=None,
                   help="exact surface-solar variable name inside the rad "
                        "files; omit to auto-detect")
    p.add_argument("--list-vars", action="store_true",
                   help="print the variables present in each archive and exit")
    p.add_argument("--chunks", type=int, nargs=3, default=list(DEFAULT_CHUNKS),
                   metavar=("TIME", "LAT", "LON"))
    p.add_argument("--max-workers", type=int, default=16)
    p.add_argument("--clevel", type=int, default=5)
    p.add_argument("--quantize-digits", type=int, default=None,
                   help="round to N decimals before compressing")
    p.add_argument("--zarr-format", type=int, choices=[2, 3], default=2,
                   help="2 = widest client compatibility (default)")
    p.add_argument("--no-deaccumulate", action="store_true",
                   help="store accumulated flux variables (e.g. rad/SSRD) "
                        "raw, without converting J m-2 -> W m-2")
    a = p.parse_args()

    names = None
    if a.vars:
        names = known if "all" in a.vars else a.vars
    meta = parse_meta(a.meta, set(known))

    if a.list_vars:
        list_vars(names, a.rad_var, root=a.root)
        return
    if not a.out:
        p.error("--out is required (or use --list-vars)")

    if a.annotate:
        if not meta:
            p.error("--annotate needs at least one --meta VAR.attr=value")
        annotate(a.out, meta)
        return

    if a.fill:
        if a.append:
            p.error("--fill and --append are different operations; use one")
        if not names:
            p.error("--fill needs an explicit --vars list")
        fill(a.out, names, root=a.root, rad_var=a.rad_var, months=a.months,
             max_workers=a.max_workers, quantize_digits=a.quantize_digits)
        return

    if a.append:
        if not names:
            p.error("--append needs an explicit --vars list (it will not "
                    "guess which variables to add)")
        append(a.out, names, root=a.root, rad_var=a.rad_var,
               replace=tuple(a.replace), meta=meta, chunks=tuple(a.chunks),
               max_workers=a.max_workers, clevel=a.clevel,
               quantize_digits=a.quantize_digits,
               deaccumulate=not a.no_deaccumulate, resume=a.resume)
        return

    if a.replace:
        p.error("--replace only applies with --append")
    if a.resume:
        p.error("--resume only applies with --append")
    if a.months:
        p.error("--months only applies with --fill")
    convert(a.out, var_names=names, rad_var=a.rad_var, root=a.root,
            chunks=tuple(a.chunks), max_workers=a.max_workers,
            clevel=a.clevel, quantize_digits=a.quantize_digits,
            zarr_format=a.zarr_format, deaccumulate=not a.no_deaccumulate,
            meta=meta)


if __name__ == "__main__":
    main()
