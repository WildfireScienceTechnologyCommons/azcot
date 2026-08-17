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
import concurrent.futures as cf
import contextlib
import datetime as dt
import os
import re

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


def build_specs(rad_var=None):
    return {
        "2T": VarSpec(
            "2T", "2T", _p_2t, src_var="2T_GDS0_SFC",
            units="K", long_name="2 metre temperature"),
        "SSRD": VarSpec(
            "SSRD", "rad", _p_rad, src_var=rad_var,
            prefer=RAD_VAR_PREFERENCE,
            long_name="surface solar irradiance",
            note=("rad files hold ~6 variables; name auto-detected or set "
                  "via --rad-var. Check units: accumulated J m-2 needs "
                  "/3600 for W m-2.")),
        "WS10": VarSpec(
            "WS10", "windclimo", _p_ws10, src_var=None,
            units="m s-1", long_name="10 metre wind speed",
            note="WS10 is precomputed speed; no u/v combination needed"),
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


@contextlib.contextmanager
def open_source(root, fs, spec, ts):
    """
    Open one hourly source file as an xarray Dataset, regardless of
    whether root is a local path or an http(s) URL.
    """
    relpath = spec.path(ts.month, ts.day, ts.hour)
    if fs is not None:
        f = fs.open(f"{root}/{relpath}", mode="rb")
        try:
            with xr.open_dataset(f) as ds:
                yield ds
        finally:
            f.close()
    else:
        with xr.open_dataset(os.path.join(root, relpath)) as ds:
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


def convert(out_path, var_names=None, rad_var=None, chunks=DEFAULT_CHUNKS,
            max_workers=16, clevel=5, quantize_digits=None, zarr_format=2,
            deaccumulate=True, verbose=True):
    """
    Build one Zarr store containing every requested variable.

    All variables must share a grid to live in one store. If the
    archives disagree this stops with a clear message rather than
    silently regridding - resampling climatology fields is a scientific
    decision, not a plumbing one.

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
    """
def convert(out_path, var_names=None, rad_var=None, root=DEFAULT_ROOT,
            chunks=DEFAULT_CHUNKS, max_workers=16, clevel=5,
            quantize_digits=None, zarr_format=2, deaccumulate=True,
            verbose=True):
    """
    Build one Zarr store containing every requested variable.

    root : local filesystem path (default: the Qumulo mount at
        /qumulo/azcot) or an http(s) URL - auto-detected.

    All variables must share a grid to live in one store. If the
    archives disagree this stops with a clear message rather than
    silently regridding - resampling climatology fields is a scientific
    decision, not a plumbing one.

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
    """
    specs = build_specs(rad_var)
    var_names = list(var_names) if var_names else list(specs)
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
        """Self-describing group label, e.g. 'grid_0p25'."""
        lats_ = info_lats[(nlat, nlon)]
        dlat = abs(float(lats_[1] - lats_[0]))
        return "grid_" + f"{dlat:.2f}".replace(".", "p")

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
    tvals = ((times - pd.Timestamp(f"{PLACEHOLDER_YEAR}-01-01"))
             // pd.Timedelta("1h")).values.astype("int64")

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

    arrays, src_vars, deaccum_by_var = {}, {}, {}
    for (nla, nlo), members in grids.items():
        grp = zroot if single_grid else zroot.create_group(group_name(nla, nlo))
        lats = info_lats[(nla, nlo)]
        lons = info[members[0]]["lons"]

        # fill_value=None on coords is essential: zarr v2 defaults it to
        # 0 and xarray then reads 0 as missing - which would silently
        # NaN out lon=0.0 (the prime meridian IS on this grid) and time
        # index 0.
        zlat = _make_array(grp, "lat", (nla,), (nla,), "float32",
                           fill_value=None)
        zlat[:] = lats
        zlat.attrs["_ARRAY_DIMENSIONS"] = ["lat"]
        zlat.attrs["units"] = "degrees_north"

        zlon = _make_array(grp, "lon", (nlo,), (nlo,), "float32",
                           fill_value=None)
        zlon[:] = lons
        zlon.attrs["_ARRAY_DIMENSIONS"] = ["lon"]
        zlon.attrs["units"] = "degrees_east"

        ztime = _make_array(grp, "time", (nt,), (nt,), "int64",
                            fill_value=None)
        ztime[:] = tvals
        ztime.attrs["_ARRAY_DIMENSIONS"] = ["time"]
        ztime.attrs["units"] = (
            f"hours since {PLACEHOLDER_YEAR}-01-01 00:00:00")
        ztime.attrs["calendar"] = "standard"
        grp.attrs["chunking_rationale"] = (
            f"chunked ({ct},{cla},{clo}) so a point timeseries costs "
            f"{int(np.ceil(nt/ct))} request(s) instead of {nt}")

        for nm in members:
            spec = specs[nm]
            with open_source(root, fs, spec, info[nm]["sample_ts"]) as ds:
                latn, lonn = find_coord_names(ds)
                sv = pick_variable(ds, spec, latn, lonn)
                src_units = ds[sv].attrs.get("units")
                src_long = ds[sv].attrs.get("long_name")
            src_vars[nm] = sv

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

            a = _make_array(grp, nm, (nt, nla, nlo), (ct, cla, clo),
                            "float32", compressor=compressor,
                            fill_value=np.nan)
            a.attrs["_ARRAY_DIMENSIONS"] = ["time", "lat", "lon"]
            a.attrs["units"] = out_units
            a.attrs["long_name"] = spec.long_name or src_long or nm
            a.attrs["source_variable"] = sv
            a.attrs["source_units"] = src_units or ""
            a.attrs["source_layout"] = f"{root}/{spec.path(1, 1, 0)}"
            a.attrs["grid"] = f"{nla}x{nlo}"
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
            arrays[nm] = a
            deaccum_by_var[nm] = deaccum_seconds

    missing = {nm: 0 for nm in var_names}
    for start in range(0, nt, ct):
        stop = min(start + ct, nt)
        block = times[start:stop]
        for nm in var_names:
            spec, sv = specs[nm], src_vars[nm]
            nla, nlo = info[nm]["nlat"], info[nm]["nlon"]
            buf = np.full((stop - start, nla, nlo), np.nan, dtype="float32")
            with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
                futs = {ex.submit(_fetch, root, fs, spec, sv, ts): i
                        for i, ts in enumerate(block)}
                for fut in cf.as_completed(futs):
                    i = futs[fut]
                    field = fut.result()
                    if field is None:
                        missing[nm] += 1
                    else:
                        buf[i] = field
            if deaccum_by_var.get(nm):
                buf = buf / deaccum_by_var[nm]
            if quantize_digits is not None:
                buf = np.round(buf, quantize_digits)
            arrays[nm][start:stop, :, :] = buf
        if verbose:
            print(f"  wrote {stop}/{nt} timesteps ({100*stop/nt:.0f}%)")

    zarr.consolidate_metadata(out_path)   # 1 metadata GET, not thousands

    if verbose:
        print(f"\nDone -> {out_path}")
        for nm in var_names:
            tag = (f" (deaccumulated /{deaccum_by_var[nm]}s -> W m-2)"
                   if deaccum_by_var.get(nm) else "")
            print(f"  {nm}: {missing[nm]} missing/unreadable hours{tag}")
        if not single_grid:
            print("\nRead with:")
            for (nla, nlo), members in grids.items():
                print(f"  xr.open_zarr(url, group='{group_name(nla, nlo)}')"
                      f"  # {', '.join(members)}")
    return out_path


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", help="output .zarr path")
    p.add_argument("--root", default=DEFAULT_ROOT,
                   help=f"source root: local path or http(s) URL "
                        f"(default: {DEFAULT_ROOT})")
    p.add_argument("--vars", nargs="*", default=None,
                   choices=["2T", "SSRD", "WS10"],
                   help="variables to include (default: all three)")
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

    if a.list_vars:
        list_vars(a.vars, a.rad_var, root=a.root)
        return
    if not a.out:
        p.error("--out is required (or use --list-vars)")

    convert(a.out, var_names=a.vars, rad_var=a.rad_var, root=a.root,
            chunks=tuple(a.chunks), max_workers=a.max_workers,
            clevel=a.clevel, quantize_digits=a.quantize_digits,
            zarr_format=a.zarr_format, deaccumulate=not a.no_deaccumulate)


if __name__ == "__main__":
    main()
