#!/usr/bin/env python3
"""
Diagnose the SSRD NaN pattern: is the missing data a fixed spatial
mask, or does it vary in time?

Motivation
----------
verify_zarr.py confirmed the Zarr store faithfully reproduces the
source netCDF (zero mismatches), but reported that ~80% of sampled
SSRD points are NaN in BOTH, while 2T and WS10 have none. The
full-timeseries sampling showed a perfectly clean per-point split
(some grid points NaN at every single hour, others valid at every
hour), which points at a static spatial mask rather than anything
time-dependent.

Worth stating plainly: this is NOT polar night. During polar night
solar irradiance is 0, not missing, and a night/season effect would
make each point partly valid rather than all-or-nothing.

This script settles it by:
  1. Computing the valid (non-NaN) mask at several widely separated
     timesteps and checking whether the mask is identical across them.
  2. Reporting the lat/lon bounding box and coverage of the valid
     region, so you can see what sub-domain actually carries data.
  3. Printing a coarse ASCII map of the mask.

Usage
-----
    python diagnose_nan_mask.py --zarr /qumulo/azcot/azcot.zarr
    python diagnose_nan_mask.py --zarr ... --var SSRD --group grid_0p10
"""

import argparse

import numpy as np
import xarray as xr
import zarr


def find_group_for(zarr_path, var):
    """Locate which group holds `var`, from the store's own attrs."""
    g = zarr.open_group(zarr_path, mode="r")
    groups_meta = g.attrs.get("groups")
    if groups_meta:
        for gname, meta in groups_meta.items():
            if var in meta["variables"]:
                return gname
    return None


def ascii_map(mask, rows=24, cols=72):
    """Coarse ASCII rendering of a 2-D boolean valid-mask."""
    nla, nlo = mask.shape
    out = []
    for r in range(rows):
        l0, l1 = int(r * nla / rows), max(int((r + 1) * nla / rows), int(r * nla / rows) + 1)
        line = []
        for c in range(cols):
            o0, o1 = int(c * nlo / cols), max(int((c + 1) * nlo / cols), int(c * nlo / cols) + 1)
            frac = mask[l0:l1, o0:o1].mean()
            line.append(" " if frac == 0 else ("." if frac < 0.5 else
                        ("+" if frac < 1.0 else "#")))
        out.append("".join(line))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--zarr", required=True)
    p.add_argument("--var", default="SSRD")
    p.add_argument("--group", default=None,
                   help="zarr group; auto-detected from store attrs if omitted")
    p.add_argument("--n-times", type=int, default=6,
                   help="how many timesteps to compare masks across")
    a = p.parse_args()

    group = a.group or find_group_for(a.zarr, a.var)
    ds = xr.open_zarr(a.zarr, group=group) if group else xr.open_zarr(a.zarr)
    da = ds[a.var]
    nt = da.sizes["time"]
    print(f"{a.var} in group '{group}': {dict(da.sizes)}")
    print(f"units: {da.attrs.get('units')}  source_variable: "
          f"{da.attrs.get('source_variable')}\n")

    # Widely separated timesteps: different months, different hours.
    idxs = sorted(set(int(round(i)) for i in
                      np.linspace(0, nt - 1, a.n_times)))

    print("=== Is the NaN mask time-invariant? ===")
    ref_mask = None
    all_same = True
    skipped = 0
    compared = 0
    for i in idxs:
        field = da.isel(time=i).values
        mask = ~np.isnan(field)
        t = str(da.time.values[i])[:16]
        frac = mask.mean()

        # An entirely-NaN timestep means that hour is absent from the
        # source archive, not that the spatial mask changed. Comparing
        # against it would be meaningless, so skip it.
        if not mask.any():
            print(f"  t[{i:5d}] {t}  valid   0.0%   (hour absent from "
                  f"archive - skipped)")
            skipped += 1
            continue

        if ref_mask is None:
            ref_mask = mask
            print(f"  t[{i:5d}] {t}  valid {frac*100:5.1f}%   (reference)")
        else:
            same = bool((mask == ref_mask).all())
            all_same &= same
            compared += 1
            ndiff = int((mask != ref_mask).sum())
            print(f"  t[{i:5d}] {t}  valid {frac*100:5.1f}%   "
                  f"identical to reference: {same}"
                  + ("" if same else f"  ({ndiff} cells differ)"))

    print()
    if ref_mask is None:
        print("  -> Every sampled timestep was entirely NaN. Increase "
              "--n-times or check the store.")
        return
    if skipped:
        print(f"  ({skipped} sampled timestep(s) skipped as absent hours)")
    if compared == 0:
        print("  -> Only one populated timestep sampled; cannot assess "
              "time-invariance. Increase --n-times.")
    elif all_same:
        print(f"  -> MASK IS TIME-INVARIANT across {compared} compared "
              f"timesteps: the same cells are missing at every one. This is a")
        print("     fixed spatial domain limit in the source rad archive, not")
        print("     a physical/diurnal/seasonal effect.")
    else:
        print("  -> Mask VARIES with time; the missing data is not a static")
        print("     domain mask. Investigate the differing timesteps above.")

    print("\n=== Where is the valid region? ===")
    lats, lons = ds.lat.values, ds.lon.values
    rows_any = ref_mask.any(axis=1)
    cols_any = ref_mask.any(axis=0)
    if not rows_any.any():
        print("  No valid data at all in the reference timestep.")
        return
    la0, la1 = np.where(rows_any)[0][[0, -1]]
    lo0, lo1 = np.where(cols_any)[0][[0, -1]]
    print(f"  valid fraction : {ref_mask.mean()*100:.2f}% of the grid")
    print(f"  lat index range: {la0}..{la1}  "
          f"-> {lats[la0]:.3f} .. {lats[la1]:.3f} degrees_north")
    print(f"  lon index range: {lo0}..{lo1}  "
          f"-> {lons[lo0]:.3f} .. {lons[lo1]:.3f} degrees_east")

    sub = ref_mask[la0:la1+1, lo0:lo1+1]
    print(f"  bounding box is {sub.mean()*100:.1f}% filled "
          f"({'a solid rectangle' if sub.all() else 'NOT solid - irregular shape'})")

    print("\n=== Valid-data map (#=all valid, +=mostly, .=some, blank=none) ===")
    print(f"  lat {lats[0]:.1f} (top) .. {lats[-1]:.1f} (bottom), "
          f"lon {lons[0]:.1f} .. {lons[-1]:.1f}")
    for line in ascii_map(ref_mask):
        print("  |" + line + "|")


if __name__ == "__main__":
    main()
