# azcot climatology tools

Reorganizes the azcot 1991-2020 hourly climatology archive into a
timeseries-optimized [Zarr](https://zarr.dev/) store, verifies the
conversion is faithful, diagnoses a data quirk found along the way, and
provides a query helper for the finished store.

## The problem this solves

The source archive is one netCDF file **per hour**, each holding a full
global grid:

```
2T/<MM>/<DD>/2T.<mon><DD><HH>.91-20clim.nc         121x1440   0.25 deg
rad/<MM>/<MM><DD><HH>.91-20radclim.nc              301x3600   0.10 deg
windclimo/<MM>/<DD>/<MM><DD><HH>.WS10.nc           121x1440   0.25 deg
```

That layout is optimal for "give me a global map at one instant" and
pessimal for "give me a timeseries at one point" - a single point-year
requires 8,784 file reads and several GB transferred to extract 8,784
numbers.

Merging into one big netCDF/HDF5 file would **not** fix this for remote
consumers: without an OPeNDAP/THREDDS server, an HTTP client has to
fetch the whole file regardless of format. Zarr solves it properly by
splitting the array into many small, independently addressable chunks,
so a client - local or remote - reads only the chunks it needs from an
ordinary static file server.

## Data quirks discovered along the way

These aren't bugs in the tooling - they're real properties of the
source archive that the tools account for explicitly rather than
silently working around:

- **Three source archives, three different directory layouts.** `2T`
  nests by day with month abbreviations in the filename; `rad` is flat
  (no day subdirectory); `windclimo` nests by day with numeric months.
  See `VarSpec` / `build_specs()` in `reorganize_to_zarr.py`.

- **Two different grids.** `2T` and `WS10` are on a 0.25° grid
  (121x1440); `SSRD` (radiation) is on a finer 0.10° grid (301x3600).
  They are **not** regridded onto a common mesh - 0.25/0.10 isn't an
  integer ratio, downsampling SSRD would throw away real resolution,
  and for point timeseries a shared grid isn't needed anyway. Each
  grid gets its own Zarr group (e.g. `grid_0p25`, `grid_0p10`),
  auto-named from actual spacing.

- **SSRD is stored as accumulated energy, not power.** The source
  variable is literally named `SSRD_GDS0_SFC_acc1h` with units
  `W m**-2 s`. Dimensionally `W*s = J`, so despite the "W" in the
  units string this is accumulated **J/m²** over the hour, not
  instantaneous irradiance. Left unconverted, values would be 3600x
  too large - and wouldn't look obviously wrong, since large numbers
  are plausible for accumulated energy. `reorganize_to_zarr.py`
  detects this from the variable name + units together (see
  `infer_deaccum_seconds()`) and converts to W/m² by default.

- **SSRD is land-masked.** ~65-80% of the SSRD grid is NaN, and the
  mask is identical at every timestep (confirmed by
  `diagnose_nan_mask.py`) - not a physical or seasonal effect. An ocean
  point returns real `2T`/`WS10` but all-NaN `SSRD` at the same
  coordinates, which looks like a broken store if you don't expect it.

- **The archive only publishes 6 of 12 months** (Jan-Mar, Oct-Dec).
  There is no year - it's a climatology - so a placeholder year (2000,
  chosen for its leap day) is used purely so datetime/pandas machinery
  works. Only month/day/hour carry real meaning.

## Scripts

### `reorganize_to_zarr.py` - build the store

Converts the three source archives into one Zarr store, one group per
grid.

```bash
./reorganize_to_zarr.py --list-vars          # inspect source variables/units first
./reorganize_to_zarr.py --out azcot.zarr      # build (local /qumulo/azcot by default)
```

Key options: `--root` (local path or http(s) URL; local skips all the
HTTP-specific workarounds and reads files natively), `--vars` (subset
of `2T`/`SSRD`/`WS10`), `--rad-var` (pin the SSRD variable name if
auto-detection picks wrong), `--chunks TIME LAT LON` (default
`2196 8 8`; larger spatial tiles trade fewer chunk-objects for bigger
per-request payloads), `--no-deaccumulate` (store SSRD raw), `--zarr-format
{2,3}` (default 2, for widest client compatibility). Works with zarr
2.18+ or 3.x via a small compatibility shim.

### `verify_zarr.py` - confirm the store matches source

Independently re-opens original netCDF files and compares actual
values - does not just trust the conversion script's own summary.

```bash
./verify_zarr.py --zarr azcot.zarr --root /qumulo/azcot \
    --n-random 300 --n-timeseries 5 --expect-no-deaccumulate
```

Three checks: random (variable, time, lat, lon) spot checks against
source; full-timeseries checks at a few fixed points (catches
time-axis misalignment, which random sampling wouldn't); and an
independent re-derivation of the SSRD deaccumulation divisor from the
source variable name (catches a class of bug that a value round-trip
alone cannot - see the script's docstring for why). Exit code 0/1, so
it can gate a build pipeline rather than just be eyeballed.

### `diagnose_nan_mask.py` - investigate a NaN pattern

Determines whether missing data in a variable is a fixed spatial mask
or something time-varying, and maps where the valid region actually
is.

```bash
./diagnose_nan_mask.py --zarr azcot.zarr --var SSRD
```

Compares the NaN mask across widely-separated timesteps (skipping
timesteps that are entirely absent from the archive, which would
otherwise look like a false "mask varies" result), reports the valid
region's bounding box, and draws a coarse ASCII map. This is how the
SSRD land mask above was confirmed rather than assumed.

### `azcot_query.py` - query the finished store

The day-to-day interface. Hides the two-grid split, the accumulated
SSRD units, and the placeholder-year time axis behind one API.

```python
from azcot_query import AzcotStore

store = AzcotStore()   # defaults to the published HTTP store; pass a
                        # local path to read a mount instead

print(store.summary())

df = store.timeseries(lat=38.9, lon=-77.0)                    # all 3 vars
df = store.timeseries(38.9, -77.0, start="01-01", end="01-31-23")
df = store.timeseries(64.8, -147.7, start="12-28", end="01-05-23")  # wraps new year

long = store.sites({"Fairbanks": (64.8, -147.7), "Reykjavik": (64.1, -21.9)})
clim = store.climatology(38.9, -77.0, by="hour")               # mean diurnal cycle

from azcot_query import plot_timeseries
plot_timeseries(df)     # month-tick x-axis; handles the Apr-Sep gap and year-wraps
```

Also usable as a CLI: `./azcot_query.py --lat 38.9 --lon -77.0 --start 01-01 --end 01-31-23`.

Notes:
- `si=True` (default) converts SSRD to W/m² using the divisor recorded
  in the store; `si=False` returns the raw accumulated value. Whichever
  you choose, the actual units returned are in `df.attrs["units"]`.
- An all-NaN result prints a note explaining the land mask rather than
  failing silently; `df.attrs["all_nan"]` lists affected variables.
- Works identically against a local path or an http(s) URL (needs
  `fsspec`+`aiohttp` for the latter).
- `plot_timeseries` auto-detects a Jupyter inline backend and detaches
  the returned figure so it isn't displayed twice.

## Requirements

```bash
pip install -U 'zarr>=2.18' numcodecs xarray numpy pandas scipy \
    h5netcdf h5py fsspec aiohttp matplotlib
```

`zarr` must be **2.18+** specifically (earlier 2.x releases use
`np.PINF`/`np.NINF`, removed in NumPy 2.0) or any 3.x. `fsspec`/`aiohttp`
are only needed when reading over HTTP; a local-path store needs
neither.

## Typical workflow

```bash
./reorganize_to_zarr.py --list-vars
./reorganize_to_zarr.py --out /srv/azcot/azcot.zarr
./verify_zarr.py --zarr /srv/azcot/azcot.zarr --root /qumulo/azcot
./diagnose_nan_mask.py --zarr /srv/azcot/azcot.zarr --var SSRD   # optional, if NaNs look surprising
```

Then, in analysis code:

```python
from azcot_query import AzcotStore
store = AzcotStore("/srv/azcot/azcot.zarr")
df = store.timeseries(lat=..., lon=...)
```
