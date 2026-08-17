#!/usr/bin/env python3
"""
Query helper for the azcot climatology Zarr store built by
reorganize_to_zarr.py.

Hides three things that would otherwise trip users up:

1. TWO GRIDS, TWO GROUPS. 2T and WS10 live on a 0.25 deg grid
   (grid_0p25); SSRD lives on a finer 0.10 deg grid (grid_0p10). They
   were deliberately not regridded onto a common mesh. This module
   opens both and picks the nearest cell in each, so you ask for a
   lat/lon once and get all three variables back - while still
   reporting the (different) actual grid cell used for each.

2. SSRD IS LAND-ONLY AND STORED ACCUMULATED.
   * Land mask: the rad source is masked to land, so an ocean point
     returns all-NaN for SSRD while 2T and WS10 return real values at
     the same coordinates. That asymmetry looks like a broken store if
     you don't expect it, so timeseries() warns explicitly.
   * Units: the store was built with --no-deaccumulate, so SSRD holds
     accumulated energy ("W m**-2 s", dimensionally J m-2 over 1h),
     NOT average irradiance. By default this module converts to W m-2
     (dividing by the accumulation period read from the array's own
     attrs). Pass si=False to get raw stored values instead. Whichever
     you choose, the units actually returned are recorded in
     df.attrs["units"].

3. THE YEAR IS MEANINGLESS. This is a 1991-2020 climatology: only
   month/day/hour carry information. The time coordinate uses a
   placeholder year (2000) purely so pandas/xarray datetime machinery
   works. Date ranges are therefore given as "MM-DD" or "MM-DD-HH",
   and a range whose end precedes its start is treated as wrapping the
   new year (e.g. "12-28" -> "01-05").

Only 6 months are published (Jan, Feb, Mar, Oct, Nov, Dec), so the
time axis contains those hours only - there are no all-NaN filler rows
for unpublished months.

Local or remote
---------------
The store path may be a local filesystem path or an http(s) URL - the
same code works for both, because Zarr is many small independently
addressable objects and a plain static file server is enough (no
OPeNDAP/THREDDS needed):

    AzcotStore()                                 # published store, over HTTP
    AzcotStore("/qumulo/azcot/azcot.zarr")       # or a local mount

HTTP access needs fsspec + aiohttp installed (pip install fsspec
aiohttp); a clear error says so if they're missing.

Consolidated metadata matters over HTTP. The store carries a single
.zmetadata document describing every array, so opening it costs ONE
request instead of one per array. This module requires it by default;
pass consolidated=False only if you know the store lacks it.

Serving the store: it is an ordinary directory tree, so any static web
server works. Two things to check:
  * The whole .zarr directory must be served, including the dotfiles
    (.zmetadata, .zgroup, .zattrs, and each array's .zarray). Some
    servers hide dot-prefixed paths by default - if .zmetadata 404s,
    that is the cause.
  * Browser-based clients (JupyterLite, zarr.js) additionally need
    CORS headers. Python/fsspec access does not.

Usage
-----
    from azcot_query import AzcotStore

    store = AzcotStore()          # defaults to the published HTTP store
    print(store.summary())

    # single point, all variables, whole archive
    df = store.timeseries(lat=38.9, lon=-77.0)

    # one month, wind only
    df = store.timeseries(38.9, -77.0, variables=["WS10"],
                          start="01-01", end="01-31-23")

    # period wrapping the new year
    df = store.timeseries(64.8, -147.7, start="12-28", end="01-05-23")

    # several sites at once -> tidy long format
    long = store.sites({"Fairbanks": (64.8, -147.7),
                        "Reykjavik": (64.1, -21.9)})

    # daily means
    df["2T"].resample("1D").mean()

Requires: xarray, zarr>=2.18, numpy, pandas
          (matplotlib only for plot_timeseries)
"""

import datetime as dt

import numpy as np
import pandas as pd
import xarray as xr
import zarr

PLACEHOLDER_YEAR = 2000

# Published store. Override with AzcotStore(path=...) or --zarr; a local
# mount (e.g. /qumulo/azcot/azcot.zarr) works identically and avoids the
# network if you are on the same filesystem.
DEFAULT_STORE = "https://scil-data.sdsc.edu/data/azcot/azcot.zarr"


def is_url(path):
    return str(path).startswith(("http://", "https://"))


def _make_mapper(url):
    """
    One fsspec mapper reused for every read, so the underlying HTTP
    session and connection pool are shared rather than rebuilt per
    group. Passing a mapper (rather than the bare URL) also works
    identically on zarr 2.x and 3.x.
    """
    try:
        import fsspec  # noqa: F401
    except ImportError:
        raise ImportError(
            "Reading a Zarr store over HTTP needs fsspec and aiohttp:\n"
            "    pip install fsspec aiohttp")
    try:
        import aiohttp  # noqa: F401
    except ImportError:
        raise ImportError(
            "fsspec needs aiohttp for http(s) access:\n"
            "    pip install aiohttp")
    import fsspec
    return fsspec.get_mapper(url)


def parse_mmddhh(s, year=PLACEHOLDER_YEAR):
    """'MM-DD' or 'MM-DD-HH' -> datetime on the placeholder year."""
    parts = str(s).split("-")
    if len(parts) == 2:
        mm, dd, hh = parts[0], parts[1], "0"
    elif len(parts) == 3:
        mm, dd, hh = parts
    else:
        raise ValueError(f"Cannot parse '{s}', expected MM-DD or MM-DD-HH")
    return dt.datetime(year, int(mm), int(dd), int(hh))


class AzcotStore:
    """Open an azcot Zarr store and query point timeseries from it."""

    def __init__(self, path=DEFAULT_STORE, si=True, consolidated=True):
        """
        path : local path OR http(s) URL to the .zarr store.
               Defaults to the published store (DEFAULT_STORE); pass a
               local path such as /qumulo/azcot/azcot.zarr to read from
               a mount instead and skip the network entirely.
        si   : convert accumulated variables to standard flux units
               (W m-2). See module docstring. Can be overridden per
               call in timeseries().
        consolidated : read the single .zmetadata document instead of
               probing every array. Leave True for HTTP - without it a
               client issues a request per array just to open the
               store, which is exactly the per-file overhead this whole
               layout exists to avoid.
        """
        self.path = str(path)
        self.default_si = si
        self.is_remote = is_url(self.path)
        self._consolidated = consolidated
        # A mapper is a MutableMapping over the store; building it once
        # keeps one HTTP session for discovery AND every group open.
        self._store = _make_mapper(self.path) if self.is_remote else self.path

        self._var_groups = self._discover_groups()
        self._groups = {}
        for g in set(self._var_groups.values()):
            self._groups[g] = self._open(g)

    def _open(self, group):
        kw = {"consolidated": self._consolidated}
        if group:
            kw["group"] = group
        try:
            return xr.open_zarr(self._store, **kw)
        except Exception as exc:
            if self._consolidated:
                raise RuntimeError(
                    f"Could not open {'group ' + group if group else 'store'} "
                    f"with consolidated metadata ({exc}). If this store was "
                    f"written without it, re-run "
                    f"zarr.consolidate_metadata(path), or construct with "
                    f"AzcotStore(..., consolidated=False) - but note that "
                    f"over HTTP the non-consolidated path is far slower."
                ) from exc
            raise

    # -- discovery -----------------------------------------------------

    def _root_attrs(self):
        """
        Read the store's root attributes.

        Done by fetching the small root metadata document directly
        rather than by opening the whole group with zarr, which saves a
        redundant round-trip per session over HTTP. Falls back to a
        real zarr open if the layout is unexpected.
        """
        import json
        for key in (".zattrs", "zarr.json"):   # v2 layout, then v3
            try:
                raw = self._store[key] if self.is_remote else None
                if raw is None:
                    import os
                    fp = os.path.join(self.path, key)
                    if not os.path.exists(fp):
                        continue
                    with open(fp, "rb") as f:
                        raw = f.read()
                doc = json.loads(raw.decode("utf-8"))
                # v3 nests user attributes under "attributes"
                return doc.get("attributes", doc)
            except (KeyError, FileNotFoundError, ValueError):
                continue
        return dict(zarr.open_group(self._store, mode="r").attrs)

    def _discover_groups(self):
        """
        Map variable -> group, read from the store's own attrs. Group
        names encode grid spacing (e.g. 'grid_0p25') and differ per
        store, so this must not be hardcoded.
        """
        attrs = self._root_attrs()
        meta = attrs.get("groups")
        if meta:
            return {v: gname for gname, m in meta.items()
                    for v in m["variables"]}
        excluded = {"lat", "lon", "time"}
        g = zarr.open_group(self._store, mode="r")
        return {v: None for v in g.array_keys() if v not in excluded}

    @property
    def variables(self):
        return sorted(self._var_groups)

    def _da(self, var):
        return self._groups[self._var_groups[var]][var]

    def _accum_seconds(self, var):
        """
        Seconds of accumulation for a variable stored raw-accumulated,
        or None if it needs no conversion (either already a flux, or
        already deaccumulated at write time).
        """
        a = self._da(var).attrs
        if a.get("deaccumulated"):
            return None          # already converted on write
        secs = a.get("deaccumulation_seconds")
        if secs:
            return int(secs)
        # Fall back to parsing the source variable name, for stores
        # written before deaccumulation metadata existed.
        sv = str(a.get("source_variable", ""))
        units = str(a.get("source_units", a.get("units", ""))).lower()
        import re
        m = re.search(r"_acc(\d+)h$", sv, re.IGNORECASE)
        looks_accumulated = (units.startswith("w m") and units.endswith("s")) \
            or units.replace("*", "").replace("-", "").startswith("j m2") \
            or units in ("j m**-2", "j m-2")
        if m and looks_accumulated:
            return int(m.group(1)) * 3600
        return None

    def summary(self):
        """Human-readable description of what's in the store."""
        lines = [f"Store: {self.path}"]
        for v in self.variables:
            da = self._da(v)
            g = self._var_groups[v]
            secs = self._accum_seconds(v)
            note = ""
            if secs:
                note = (f"  [stored accumulated over {secs//3600}h; "
                        f"si=True divides by {secs} -> W m-2]")
            lines.append(
                f"  {v:6s} group={g or '<root>':10s} "
                f"grid={da.sizes['lat']}x{da.sizes['lon']} "
                f"units={da.attrs.get('units','?')}{note}")
        t = self._groups[self._var_groups[self.variables[0]]].time
        months = sorted(pd.DatetimeIndex(t.values).month.unique())
        lines.append(f"  time: {len(t)} hourly steps, months {months} "
                     f"(year is a placeholder and carries no meaning)")
        return "\n".join(lines)

    # -- querying ------------------------------------------------------

    def _time_mask(self, times, start, end):
        """Boolean mask, plus whether the range wraps the year boundary."""
        if start is None and end is None:
            return np.ones(len(times), dtype=bool), False
        s = parse_mmddhh(start) if start else None
        e = parse_mmddhh(end) if end else None
        if s is not None and e is not None and e < s:
            # Wrap across the new year: keep >= start OR <= end.
            return ((times >= s) | (times <= e)), True
        mask = np.ones(len(times), dtype=bool)
        if s is not None:
            mask &= times >= s
        if e is not None:
            mask &= times <= e
        return mask, False

    def timeseries(self, lat, lon, variables=None, start=None, end=None,
                   si=None, warn=True):
        """
        Point timeseries at the cell nearest (lat, lon).

        lat, lon  : degrees; lon in -180..180 (NOT 0..360)
        variables : subset to fetch (default: all)
        start,end : 'MM-DD' or 'MM-DD-HH'; end < start wraps the year
        si        : override the store default for unit conversion
        warn      : print a note if a variable is entirely NaN (the
                    usual cause is the SSRD land mask - see below)

        Returns a DataFrame indexed by time, one column per variable.
        df.attrs holds 'units', 'grid_cell' (the actual lat/lon used
        per variable - these differ between the two grids) and
        'all_nan' (variables that returned no data at all).
        """
        si = self.default_si if si is None else si
        variables = list(variables) if variables else self.variables

        cols, units, cells, all_nan = {}, {}, {}, []
        index = None

        for v in variables:
            if v not in self._var_groups:
                raise KeyError(f"{v!r} not in store; have {self.variables}")
            da = self._da(v)
            sel = da.sel(lat=lat, lon=lon, method="nearest")
            times = pd.DatetimeIndex(sel.time.values)
            mask, wrapped = self._time_mask(times, start, end)
            vals = np.asarray(sel.values, dtype="float64")[mask]
            t = times[mask]

            u = da.attrs.get("units", "")
            secs = self._accum_seconds(v)
            if si and secs:
                vals = vals / secs
                u = "W m-2"

            if index is None:
                index = t
            elif not index.equals(t):
                # Both groups are written with an identical time axis,
                # so this should not happen; fail loudly if it ever does
                # rather than silently misaligning columns.
                raise RuntimeError(
                    f"time axis mismatch between variables: {v} differs")

            cols[v] = vals
            units[v] = u
            cells[v] = (float(sel.lat), float(sel.lon))
            if vals.size and np.all(np.isnan(vals)):
                all_nan.append(v)

        df = pd.DataFrame(cols, index=index)
        df.index.name = "time"
        df.attrs["units"] = units
        df.attrs["grid_cell"] = cells
        df.attrs["all_nan"] = all_nan
        df.attrs["requested"] = (lat, lon)
        # A wrap query (e.g. 12-28 -> 01-05) selects rows at BOTH ends
        # of the placeholder year. The frame stays sorted so pandas
        # time ops (resample, rolling) keep working on a monotonic
        # index; plotting reorders for display only, using this flag
        # plus the start month.
        df.attrs["wrapped"] = wrapped
        df.attrs["wrap_start_month"] = (parse_mmddhh(start).month
                                        if wrapped else None)

        if warn and all_nan:
            for v in all_nan:
                extra = ""
                if v == "SSRD":
                    extra = (" SSRD is masked to land in the source archive, "
                             "so ocean points have no data. 2T/WS10 are not "
                             "masked and will still return values here.")
                print(f"NOTE: {v} is entirely NaN at "
                      f"({lat}, {lon}) -> nearest cell {cells[v]}.{extra}")
        return df

    def sites(self, locations, variables=None, start=None, end=None,
              si=None, warn=True):
        """
        Same as timeseries() for several named locations, returned in
        tidy long format: columns [site, lat, lon, time, variable,
        value]. Handy for groupby/seaborn work.

        locations : {name: (lat, lon)} or [(lat, lon), ...]
        """
        if not isinstance(locations, dict):
            locations = {f"site_{i}": tuple(p) for i, p in enumerate(locations)}

        frames = []
        units = {}
        for name, (la, lo) in locations.items():
            df = self.timeseries(la, lo, variables, start, end, si, warn)
            units.update(df.attrs["units"])
            m = df.reset_index().melt(id_vars="time", var_name="variable",
                                      value_name="value")
            m.insert(0, "site", name)
            m["lat"] = la
            m["lon"] = lo
            frames.append(m)

        out = pd.concat(frames, ignore_index=True)
        out.attrs["units"] = units
        return out

    def climatology(self, lat, lon, variables=None, by="dayofyear", si=None):
        """
        Aggregate a point timeseries to a mean annual cycle.

        by : 'dayofyear' (daily means), 'hour' (mean diurnal cycle),
             or 'month'.
        """
        df = self.timeseries(lat, lon, variables, si=si, warn=False)
        idx = df.index
        key = {"dayofyear": idx.dayofyear, "hour": idx.hour,
               "month": idx.month}[by]
        out = df.groupby(key).mean()
        out.index.name = by
        out.attrs.update(df.attrs)
        return out


def _month_starts(index):
    """Positions in `index` where a new month begins (first row of each)."""
    months = pd.Series(index.month, index=range(len(index)))
    change = months.ne(months.shift())
    return list(change[change].index)


def _gap_positions(index, gap=pd.Timedelta("2h")):
    """Positions where the time index jumps (missing months, or a wrap)."""
    if len(index) < 2:
        return []
    d = pd.Series(index).diff()
    return [i for i in range(1, len(index)) if d.iloc[i] > gap]


def _apply_month_ticks(ax, index, compress):
    """
    Month ticks that stay readable for this data.

    Two axis modes, because the published archive is discontinuous
    (only Jan-Mar and Oct-Dec exist) and a wrap-around query returns
    January and December sitting in the same placeholder year:

    compress=True  - x is row position, so gaps take no width. Ticks
                     sit at the first row of each month. Discontinuities
                     are drawn as dashed vertical rules so the jump is
                     visible rather than hidden.
    compress=False - real datetime axis; true spacing is preserved and
                     missing months show as blank stretches.
    """
    import matplotlib.dates as mdates

    span_days = (index[-1] - index[0]) / pd.Timedelta("1D") if len(index) > 1 else 0

    if compress:
        starts = _month_starts(index)
        ax.set_xticks(starts)
        ax.set_xticklabels([index[i].strftime("%b") for i in starts])
        # Short span: month labels alone are useless, add the day.
        if span_days <= 62 and len(starts) <= 2:
            ticks = list(range(0, len(index), max(1, len(index) // 8)))
            ax.set_xticks(ticks)
            ax.set_xticklabels([index[i].strftime("%b %-d") for i in ticks])
        for g in _gap_positions(index):
            ax.axvline(g - 0.5, color="0.5", linestyle="--", linewidth=0.8)
        ax.set_xlim(0, len(index) - 1)
        return

    if span_days > 45:
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b"))
        ax.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=(1, 15)))
    elif span_days > 5:
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=max(1, int(span_days // 8))))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d"))
    else:
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d %H:%M"))


def _is_inline_backend():
    """
    True when matplotlib will auto-display figures at end of cell
    (Jupyter's inline/ipympl backends), as opposed to a script backend
    where the user calls plt.show() themselves.
    """
    import matplotlib
    b = matplotlib.get_backend().lower()
    return ("inline" in b) or ("ipympl" in b) or ("widget" in b) or ("nbagg" in b)


def plot_timeseries(df, variables=None, title=None, compress_gaps="auto",
                    show=None):
    """
    Multi-panel plot of a timeseries() result, with month ticks.

    compress_gaps : True, False, or "auto" (default). The archive
        publishes only 6 months, so a real datetime axis leaves roughly
        half the plot empty; "auto" switches to a gap-free position
        axis when more than 10% of the plotted span is missing, and
        marks each discontinuity with a dashed rule. Set False to force
        a true-to-scale datetime axis.

    show : controls the notebook double-render problem. In Jupyter the
        inline backend auto-displays every figure created during a
        cell, AND the returned Figure is displayed again by its repr -
        so a bare `plot_timeseries(df)` renders twice. Default (None)
        detects an inline backend and detaches the figure from pyplot's
        registry, leaving exactly one render (from the return value)
        while still handing back a usable Figure for savefig(). Pass
        show=True to force pyplot to keep it (script use with
        plt.show()), or show=False to always detach.

    Returns the Figure.
    """
    import matplotlib.pyplot as plt

    variables = variables or [c for c in df.columns
                              if c not in df.attrs.get("all_nan", [])]
    if not variables:
        raise ValueError("nothing to plot - all requested variables are NaN")
    if len(df) == 0:
        raise ValueError("nothing to plot - the selection is empty")

    index = df.index
    # A wrap query returns e.g. Jan 1-5 AND Dec 28-31 in one placeholder
    # year, sorted Jan-first. Plotted as-is that draws the period
    # backwards, so roll it to lead with the pre-new-year part. Display
    # only - the caller's frame is untouched.
    if df.attrs.get("wrapped") and df.attrs.get("wrap_start_month"):
        lead = index.month >= df.attrs["wrap_start_month"]
        if lead.any() and (~lead).any():
            df = pd.concat([df[lead], df[~lead]])
            index = df.index

    if compress_gaps == "auto":
        if df.attrs.get("wrapped"):
            # Reordered wrap: the index deliberately runs Dec -> Jan, so
            # a real datetime axis would span backwards. Position axis.
            compress = True
        elif len(index) > 1:
            span = (index[-1] - index[0]) / pd.Timedelta("1h") + 1
            compress = span > 0 and (span - len(index)) / span > 0.10
        else:
            compress = False
    else:
        compress = bool(compress_gaps)

    x = range(len(index)) if compress else index

    fig, axes = plt.subplots(len(variables), 1, sharex=True,
                             figsize=(11, 2.6 * len(variables)))
    if len(variables) == 1:
        axes = [axes]
    for ax, v in zip(axes, variables):
        ax.plot(x, df[v].values, linewidth=0.7)
        ax.set_ylabel(f"{v}\n({df.attrs['units'].get(v,'')})")
        ax.grid(alpha=0.3)

    _apply_month_ticks(axes[-1], index, compress)

    lat, lon = df.attrs.get("requested", (None, None))
    axes[0].set_title(title or f"azcot climatology at {lat}, {lon}")
    xlabel = "month (climatology; year is a placeholder)"
    if compress and _gap_positions(index):
        xlabel += "  - dashed rules mark gaps where months are unpublished"
    axes[-1].set_xlabel(xlabel)
    fig.tight_layout()

    # Detach from pyplot's registry so the inline backend does not
    # ALSO auto-display it at end of cell. The Figure stays fully
    # usable (savefig, further axes tweaks) and renders once via the
    # return value's repr.
    keep = _is_inline_backend() is False if show is None else bool(show)
    if not keep:
        plt.close(fig)
    return fig


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description="Quick CLI over the azcot store")
    p.add_argument("--zarr", default=DEFAULT_STORE,
                   help=f"store path or URL (default: {DEFAULT_STORE})")
    p.add_argument("--lat", type=float)
    p.add_argument("--lon", type=float)
    p.add_argument("--vars", nargs="*", default=None)
    p.add_argument("--start", default=None)
    p.add_argument("--end", default=None)
    p.add_argument("--raw", action="store_true",
                   help="return stored values without converting "
                        "accumulated variables to W m-2")
    p.add_argument("--out", default=None, help="write CSV here")
    a = p.parse_args()

    store = AzcotStore(a.zarr, si=not a.raw)
    print(store.summary())
    if a.lat is None or a.lon is None:
        raise SystemExit(0)

    df = store.timeseries(a.lat, a.lon, a.vars, a.start, a.end)
    print()
    print("grid cell used per variable:", df.attrs["grid_cell"])
    print("units:", df.attrs["units"])
    print(df.describe())
    if a.out:
        df.to_csv(a.out)
        print(f"\nwrote {len(df)} rows to {a.out}")
