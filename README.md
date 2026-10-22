# azcot climatology tools

Reorganizes the azcot 1991-2020 hourly climatology archive into a
timeseries-optimized [Zarr](https://zarr.dev/) store, verifies the
conversion is faithful, and provides a query helper for the finished
store.

## The problem this solves

The source archive is one netCDF file **per hour**, each holding a full
global grid, in several different layouts:

```
2T/<MM>/<DD>/2T.<mon><DD><HH>.91-20clim.nc          121x1440   0.25 deg
rad/<MM>/<MM><DD><HH>.91-20radclim.nc               301x3600   0.10 deg
windclimo/<MM>/<DD>/<MM><DD><HH>.WS10.nc            121x1440   0.25 deg
STL1/<MM>/STL1.<mon><DD><HH>.91-20clim.nc           121x1440   0.25 deg
sd/<MM>/SD.<mon><DD><HH>.91-20clim.nc               121x1440   0.25 deg
skt/<MM>/SKT.<mon><DD><HH>.91-20clim.nc             121x1440   0.25 deg
wct/<MM>/<MM><DD><HH>.91-20wctclim.nc               121x1440   0.25 deg
```

That layout is optimal for "give me a global map at one instant" and
pessimal for "give me a timeseries at one point": a single point-year
requires 8,784 file reads and several GB transferred to extract 8,784
numbers.

Merging into one big netCDF/HDF5 file would **not** fix this for remote
consumers: without an OPeNDAP/THREDDS server, an HTTP client has to
fetch the whole file regardless of format. Zarr solves it properly by
splitting the array into many small, independently addressable chunks,
so a client - local or remote - reads only the chunks it needs from an
ordinary static file server.

## Variables in the store

| variable | meaning | units | group (grid) | months with data |
|---|---|---|---|---|
| `2T` | 2 m temperature | K | `grid_0p25` | 1 2 3 10 11 12 |
| `WS10` | 10 m wind speed | m s-1 | `grid_0p25` | 1 2 3 10 11 12 |
| `SSRD` | surface solar radiation, downwards | see below | `grid_0p10` | 1 2 3 10 11 12 |
| `STL1` | soil temperature level 1 | K | `grid_0p25` | 1 2 3 10 12 (no Nov) |
| `SD` | snow depth | m of water equivalent | `grid_0p25` | 1 2 3 11 12 (no Oct) |
| `SKT` | skin temperature | K | `grid_0p25` | 1 2 3 11 12 (no Oct) |
| `WCT` | unverified | **unverified** | `grid_0p25` | 1 2 3 10 11 12 |

The months `STL1`, `SD` and `SKT` lack are months whose tarballs are not
in the delivered archive, not a conversion fault. They are stored as NaN
and recorded in each array's `months_with_data` attribute.

## Data quirks

These are real properties of the source archive that the tools account
for explicitly rather than silently working around:

- **Several different directory layouts.** `2T` and `windclimo` nest by
  day; `rad`, `STL1`, `sd`, `skt` and `wct` are flat within a month
  directory; month appears as an abbreviation in some filenames and as a
  number in others. See `VarSpec` / `build_specs()`.

- **Two grids.** Everything is 0.25 deg (121x1440) except `SSRD`, which
  is 0.10 deg (301x3600). They are **not** regridded onto a common mesh:
  0.25/0.10 isn't an integer ratio, downsampling SSRD would discard real
  resolution, and for point timeseries a shared grid isn't needed. Each
  grid gets its own Zarr group, auto-named from its spacing. A new
  variable is placed in an existing group only if its lat/lon match that
  group's exactly.

- **SSRD is accumulated energy, not power.** The source variable is
  `SSRD_GDS0_SFC_acc1h` with units `W m**-2 s`. Dimensionally `W*s = J`,
  so despite the "W" this is accumulated J/m2 over the hour. Unconverted
  it is 3600x too large - and doesn't look wrong, since big numbers are
  plausible for accumulated energy. A build converts it to W m-2 by
  default (`infer_deaccum_seconds()` requires both the `_acc1h` name tag
  and energy-dimensioned units). A store built with `--no-deaccumulate`
  holds the raw values; the array's `deaccumulated` / `units` attributes
  say which you have, and `azcot_query` converts on read when needed.

- **SSRD is land-masked.** About 65% of the SSRD grid is NaN, identically
  at every timestep (confirmed with `diagnose_nan_mask.py`). An ocean
  point has real `2T`/`WS10` but all-NaN `SSRD`.

- **`WCT` has no units and no long_name** in its source files, and a NaN
  fill value rather than the 1e20 the other families use. Do not assume a
  unit. `check_wct_units.py` tests it against 2T and WS10; until the
  units are established the array records `units_note` instead of a
  value.

- **`STL1`'s variable is `STL1_GDS0_DBLY`**, not the `..._SFC` pattern of
  the others. Variable names are auto-detected rather than assumed.

- **Placeholder year.** There is no year - it's a climatology - so 2000
  (chosen for its leap day) is used purely so datetime machinery works.
  Only month, day and hour mean anything.

- **Raw per-year files are not converted.** The `2t`, `wind` and `rad`
  tarballs are ~300 GB and consist mostly of individual-year files (one
  file per year per hour), the inputs the climatologies were averaged
  from. That is a different product - it has a real time axis - and does
  not belong on this placeholder-year axis. Their climatology is already
  in the store. Using the raw files (e.g. for interannual variability)
  would need a separate store designed for 30x more time steps.

## Scripts

### `inspect_archive.py` - look before converting

Every archive added so far has differed from the last in ways filenames
didn't reveal. This inventories a directory of tarballs and reports, per
family: months present, file layout (climatology vs raw files show up as
separate patterns), and for one sample file read **in memory** (nothing
is extracted): grid, variables, units, fill value, valid fraction, and
whether anything looks like the accumulated-flux trap. With `--zarr` it
also reports whether each grid matches an existing group.

```bash
./inspect_archive.py archive/ --zarr /qumulo/azcot/azcot.zarr --out report.txt
```

Shell scripts in the directory are not printed (download scripts often
hold credentials); credential-looking lines in text files are hidden.

### `reorganize_to_zarr.py` - build, append, annotate

**Build** a store from scratch (this overwrites the output path):

```bash
./reorganize_to_zarr.py --list-vars
./reorganize_to_zarr.py --out azcot.zarr                 # 2T SSRD WS10
./reorganize_to_zarr.py --out azcot.zarr --vars all      # everything
```

**Append** variables to a store that already exists, leaving what is in
it untouched:

```bash
./reorganize_to_zarr.py --out azcot.zarr --append --vars STL1 SD SKT WCT --max-workers 4
```

Append refuses, and changes nothing, if a variable is already present
(`--replace NAME` overwrites deliberately), if the source has a month the
store's time axis lacks, or if the store isn't one this tool made. A new
variable joins the existing group its grid matches; otherwise it gets a
new group, with a disambiguated name if the spacing collides with an
existing group of a different shape. The registry attribute and
consolidated metadata are updated only after every new variable is
written, so an interrupted run leaves invisible arrays, not a store that
advertises data it lacks.

**After an interruption.** If an append of *new* variables dies (crash,
Ctrl-C, lost mount), the store is unaffected: they are registered only
once every one has been written, so what it had written is invisible.
Re-run the same command; leftovers of the variables being added are
deleted and rewritten. (`--replace` is different: the variable is
un-registered first, then its data is deleted and rewritten, so an
interruption leaves it *absent* rather than broken, and it is unavailable
to readers meanwhile.) Add `--resume` to keep the variables that finished completely
and redo only the unfinished one. A leftover is kept only if it carries
its completion marker and matches the request (same source variable,
shape and deaccumulation option); otherwise it is redone.

**Fill** months that arrived after a variable was added (a tarball for a
month the store had as NaN). It writes only those months, in place,
deleting nothing, and checks each against what is already stored (units
and grid must match) before writing:

```bash
tar -xzf archive/novSTL1.tar.gz -C /qumulo/azcot/STL1       # extract first
./reorganize_to_zarr.py --out azcot.zarr --fill --vars STL1 SD SKT --max-workers 4
./verify_zarr.py --zarr azcot.zarr --root /qumulo/azcot --vars STL1 SD SKT --expect-no-deaccumulate
```

Fill *before* verifying: until then the store holds NaN where source
files now exist, which the verifier correctly reports as a mismatch. Use
`--fill` rather than `--replace` for late months: it reads only the new
months and can't leave a variable missing. Arrays built before coverage
was recorded (the original `2T`, `WS10`, `SSRD`) need `--months`.

**Annotate** attributes without writing any data:

```bash
./reorganize_to_zarr.py --out azcot.zarr --annotate --meta WCT.units=degF
```

Other options: `--root` (local path or http(s) URL), `--chunks TIME LAT
LON` (default `2196 8 8`), `--no-deaccumulate`, `--rad-var`, `--meta
VAR.attr=value` (also usable when building or appending).

**Extracting the new tarballs.** Each family's members are already laid
out as `<MM>/<file>`, so extract each into its own directory under the
source root:

```bash
cd /qumulo/azcot
for fam in STL1 sd skt wct; do
  mkdir -p "$fam"
  for t in archive/*"$fam".tar.gz archive/*"$fam".tgz; do
    [ -e "$t" ] && tar -xzf "$t" -C "$fam"
  done
done
```

### `verify_zarr.py` - confirm the store matches source

Independently re-opens the original netCDF files and compares actual
values; it does not trust the converter's own summary.

```bash
./verify_zarr.py --zarr azcot.zarr --root /qumulo/azcot \
    --n-random 300 --n-timeseries 5 --expect-no-deaccumulate
```

Random spot checks, full-timeseries checks at fixed points (catches
time-axis misalignment), and an independent re-derivation of the SSRD
deaccumulation divisor. Exit code 0/1.

### `check_wct_units.py` - establish the WCT units from the data

Assumes WCT is a wind chill temperature (an inference from the name) and
combines three kinds of evidence against 2T and WS10 at the same hours:
a formula-free constraint (a wind chill cannot exceed the air temperature
in the right unit, which rules units out); a fit against two known
definitions, the modern (2001) one and the pre-2001 NWS (Siple-Passel)
one, each converted to each unit, where only the right (formula, unit)
pair gives slope ~1 and intercept ~0; and an identity test (many
implementations set WCT = air temperature where wind chill is undefined,
which names the unit with no formula). The pass bar is fixed:
`|slope-1| <= 0.10`, intercept within 3 C-equivalent, `R^2 >= 0.80`. It
reports **inconclusive** rather than guessing when the evidence does not
agree on exactly one unit. If you record a unit it infers, say so:
`--meta WCT.units=degF --meta "WCT.units_source=inferred from data; not
confirmed by the producer"`.

```bash
./check_wct_units.py --root /qumulo/azcot --zarr /qumulo/azcot/azcot.zarr
```

### `diagnose_nan_mask.py` - investigate a NaN pattern

Determines whether missing data in a variable is a fixed spatial mask or
time-varying, and maps where the valid region is.

```bash
./diagnose_nan_mask.py --zarr azcot.zarr --var SSRD
```

### `azcot_query.py` - query the finished store

The day-to-day interface: hides the grid split, the SSRD units, the
placeholder year and the per-variable month gaps.

```python
from azcot_query import AzcotStore

store = AzcotStore()                 # published HTTP store; or pass a path
print(store.summary())               # units, grids, months each variable lacks

df = store.timeseries(lat=38.9, lon=-77.0)
df = store.timeseries(64.8, -147.7, start="12-28", end="01-05-23")  # wraps new year
long = store.sites({"Fairbanks": (64.8, -147.7), "Reykjavik": (64.1, -21.9)})
clim = store.climatology(38.9, -77.0, by="hour")

from azcot_query import plot_timeseries
plot_timeseries(df)                  # month ticks; handles the gaps and year-wraps
```

`df.attrs` records the units returned, the actual grid cell used per
variable (the two grids differ), `all_nan` variables, and
`months_missing` - months inside your window for which a variable has no
source data. A note is printed when that happens. Stores built before
coverage was recorded report coverage as unknown rather than claiming it
is complete.

## Troubleshooting

**`Segmentation fault (core dumped)` or `free(): invalid size` during a
build or append.** Seen with a variable delivered as NetCDF-4 (HDF5)
files. Opening those from several threads at once through the netCDF4
engine corrupts memory in netCDF-C/HDF5; it is not a Python error and
leaves no message. Classic NetCDF-3 files (what the NCL-written families
use) are unaffected. The converter now opens HDF5-format files under a
lock, so `--max-workers` helps classic sources but not HDF5 ones. Check
which kind a file is:

```bash
python3 -c "print(open('/qumulo/azcot/wct/01/010100.91-20wctclim.nc','rb').read(8))"
# b'\x89HDF\r\n\x1a\n'  -> HDF5 / NetCDF-4          b'CDF\x01' / b'CDF\x02' -> classic
```

If a native crash still happens, the converter prints every thread's
Python stack at the moment of the crash (`faulthandler`); that shows
which function and file it died in. Then re-run the same command, with
`--resume` to skip variables that had finished.

**`ContainsArrayError` on a repeat append** meant leftovers from an
interrupted run blocked the re-run. That is now handled; use a current
`reorganize_to_zarr.py`.

**Grouped stores only.** `--append`, `--fill` and `check_wct_units.py`
need a store with more than one grid group, as the real one has. A store
built from a single grid is written flat at the root, without groups, and
they refuse it.

## Requirements

```bash
pip install -U 'zarr>=2.18' numcodecs xarray numpy pandas scipy \
    h5netcdf h5py fsspec aiohttp matplotlib
```

- `zarr` must be **2.18+** (earlier 2.x uses `np.PINF`/`np.NINF`, removed
  in NumPy 2.0) or any 3.x.
- **zarr 2.x needs an older xarray.** Tested: `xarray 2026.7.0` with
  `zarr 2.18.3` reads these stores; `xarray 2026.9.0` with `zarr 2.18.3`
  does not (`open_group() got an unexpected keyword argument
  'zarr_format'`, even on a store that was never modified). If you stay
  on zarr 2, pin xarray (`pip install 'xarray<=2026.7.0'`); otherwise
  move to zarr 3. The tools turn this error into a readable message.
- `fsspec`/`aiohttp` are only needed when reading over HTTP.

## Typical workflow

First build:

```bash
./reorganize_to_zarr.py --list-vars
./reorganize_to_zarr.py --out /srv/azcot/azcot.zarr
./verify_zarr.py --zarr /srv/azcot/azcot.zarr --root /qumulo/azcot
```

Adding more data later:

```bash
./inspect_archive.py archive/ --zarr /srv/azcot/azcot.zarr --out report.txt
# extract the climatology families (see above), then:
./reorganize_to_zarr.py --out /srv/azcot/azcot.zarr --append --vars STL1 SD SKT WCT
./verify_zarr.py --zarr /srv/azcot/azcot.zarr --root /qumulo/azcot --expect-no-deaccumulate
./check_wct_units.py --root /qumulo/azcot --zarr /srv/azcot/azcot.zarr
./reorganize_to_zarr.py --out /srv/azcot/azcot.zarr --annotate --meta WCT.units=<unit>
```
