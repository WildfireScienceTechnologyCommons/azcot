#!/usr/bin/env python3
"""
Verify that a Zarr store produced by reorganize_to_zarr.py contains
exactly the same data as the original per-hour netCDF archive it was
built from.

This does NOT trust the conversion script's own "0 missing hours"
summary - it independently re-opens source files and compares actual
values. Two checks are run:

  1. RANDOM SPOT CHECKS: many random (variable, time, lat, lon) points,
     each compared against the real source file for that hour. Catches
     spatial/indexing bugs (wrong lat/lon index, transposed axes) and
     value bugs (wrong scaling, deaccumulation applied incorrectly).

  2. FULL TIMESERIES CHECKS: for a few fixed grid points, every single
     hour in the zarr store is compared against source (or checked to
     be NaN if the source file for that hour is genuinely absent).
     Catches time-axis misalignment (an off-by-one in month/day/hour
     decoding would show up as a systematic shift, not random noise).

Deaccumulation-aware: if reorganize_to_zarr.py applied the accumulated
-> average-flux conversion (see infer_deaccum_seconds in that script),
the array carries deaccumulated=True and deaccumulation_seconds in its
attrs. This script reads those attrs and reverses the conversion
before comparing, so it works correctly whether or not
--no-deaccumulate was used for the run being checked.

Path construction reuses build_specs()/VarSpec/open_source() from
reorganize_to_zarr.py rather than re-deriving the URL patterns, since
those are pure, side-effect-free functions already exercised by a
successful --list-vars run against the real archive; duplicating them
here would risk the copies silently drifting apart, which is a worse
failure mode than sharing them.

Usage
-----
    python verify_zarr.py --zarr /qumulo/azcot/azcot.zarr \\
        --root /qumulo/azcot --n-random 300 --n-timeseries 5

Exits 0 if everything matches, 1 if any mismatch is found (so it can
be used as a CI/cron gate, not just eyeballed).
"""

import argparse
import random
import sys

import numpy as np
import xarray as xr
import zarr

from reorganize_to_zarr import (
    DEFAULT_ROOT, PLACEHOLDER_YEAR, build_specs, find_coord_names, make_fs,
    open_source,
)


def discover_var_groups(zarr_path):
    """
    Map each variable name to the zarr group it lives in, read from the
    store's own attrs (written by reorganize_to_zarr.py's convert()) -
    group names are self-describing based on actual grid spacing
    (e.g. "grid_0p25") and differ per store, so this must not be
    hardcoded. Returns {varname: group_name_or_None}; None means the
    variable sits flat at the store root (the single-grid case).
    """
    g = zarr.open_group(zarr_path, mode="r")
    groups_meta = g.attrs.get("groups")
    if groups_meta:
        mapping = {}
        for gname, meta in groups_meta.items():
            for v in meta["variables"]:
                mapping[v] = gname
        return mapping
    excluded = {"lat", "lon", "time"}
    return {v: None for v in g.array_keys() if v not in excluded}


def _open_group(zarr_path, group):
    return xr.open_zarr(zarr_path, group=group) if group else xr.open_zarr(zarr_path)


def _decode_time(t64):
    """numpy datetime64 (placeholder year) -> (month, day, hour)."""
    ts = np.datetime64(t64, "h").astype(object)  # -> python datetime
    return ts.month, ts.day, ts.hour


def _source_value(root, fs, spec, mm, dd, hh, lat_idx, lon_idx, src_var):
    """
    Read one raw value straight from the original source file, or
    return "missing" (a sentinel, distinct from None-as-error) if that
    hour's file genuinely doesn't exist. Returns (value_or_None, existed).
    """
    import datetime as dt
    ts = dt.datetime(PLACEHOLDER_YEAR, mm, dd, hh)
    try:
        with open_source(root, fs, spec, ts) as ds:
            latn, lonn = find_coord_names(ds)
            v = float(ds[src_var].isel({latn: lat_idx, lonn: lon_idx}).values)
            return v, True
    except Exception:
        return None, False


def _reverse_deaccum(zarr_value, attrs):
    """If this array was deaccumulated on write, undo it for comparison."""
    if attrs.get("deaccumulated"):
        secs = attrs.get("deaccumulation_seconds")
        return zarr_value * secs
    return zarr_value


def random_spot_checks(zarr_path, root, fs, specs, n_random, atol, rtol,
                       seed, verbose):
    rng = random.Random(seed)
    results = {"checked": 0, "mismatches": [], "nan_ok": 0, "value_ok": 0,
               "nan_value_ok": 0, "per_var": {}}

    var_groups = discover_var_groups(zarr_path)
    opened = {g: _open_group(zarr_path, g) for g in set(var_groups.values())}

    for nm, spec in specs.items():
        if nm not in var_groups:
            if verbose:
                print(f"  {nm}: not found in zarr store, skipping")
            continue
        grp = var_groups[nm]
        ds = opened[grp]
        da = ds[nm]
        nt, nla, nlo = da.sizes["time"], da.sizes["lat"], da.sizes["lon"]
        sv = da.attrs.get("source_variable")
        attrs = da.attrs
        stats = results["per_var"].setdefault(
            nm, {"n": 0, "nan_both": 0, "real": 0})

        for _ in range(n_random):
            ti = rng.randrange(nt)
            li = rng.randrange(nla)
            oi = rng.randrange(nlo)

            zval = float(da.isel(time=ti, lat=li, lon=oi).values)
            mm, dd, hh = _decode_time(da.time.values[ti])

            sval, existed = _source_value(root, fs, spec, mm, dd, hh, li, oi, sv)
            results["checked"] += 1
            stats["n"] += 1

            if not existed:
                # The hour's file is genuinely absent from the archive.
                if np.isnan(zval):
                    results["nan_ok"] += 1
                else:
                    results["mismatches"].append(
                        (nm, mm, dd, hh, li, oi,
                         f"source file missing but zarr has non-NaN value {zval}"))
                continue

            # The file exists - but the VALUE inside it may still be NaN,
            # decoded from the source _FillValue (1e20). NaN on both
            # sides is a correct match, not a mismatch.
            src_nan = sval is None or np.isnan(sval)
            z_nan = np.isnan(zval)

            if src_nan and z_nan:
                results["nan_value_ok"] += 1
                stats["nan_both"] += 1
                continue
            if z_nan and not src_nan:
                results["mismatches"].append(
                    (nm, mm, dd, hh, li, oi,
                     f"zarr is NaN but source has real value {sval}"))
                continue
            if src_nan and not z_nan:
                results["mismatches"].append(
                    (nm, mm, dd, hh, li, oi,
                     f"source is NaN (fill) but zarr has value {zval}"))
                continue

            stats["real"] += 1
            expected = _reverse_deaccum(zval, attrs)
            if np.isclose(expected, sval, atol=atol, rtol=rtol):
                results["value_ok"] += 1
            else:
                results["mismatches"].append(
                    (nm, mm, dd, hh, li, oi,
                     f"zarr={zval} (reversed={expected}) vs source={sval}, "
                     f"diff={abs(expected-sval):.6g}"))

        if verbose:
            print(f"  {nm}: sampled {n_random} points from {grp}")

    return results


def full_timeseries_checks(zarr_path, root, fs, specs, n_points, atol, rtol,
                           seed, verbose):
    rng = random.Random(seed + 1)
    results = {"checked": 0, "mismatches": [], "nan_ok": 0, "value_ok": 0,
               "nan_value_ok": 0, "per_var": {}}

    var_groups = discover_var_groups(zarr_path)
    opened = {g: _open_group(zarr_path, g) for g in set(var_groups.values())}

    for nm, spec in specs.items():
        if nm not in var_groups:
            continue
        grp = var_groups[nm]
        ds = opened[grp]
        da = ds[nm]
        nt, nla, nlo = da.sizes["time"], da.sizes["lat"], da.sizes["lon"]
        sv = da.attrs.get("source_variable")
        attrs = da.attrs
        stats = results["per_var"].setdefault(
            nm, {"n": 0, "nan_both": 0, "real": 0})

        for _ in range(n_points):
            li = rng.randrange(nla)
            oi = rng.randrange(nlo)
            if verbose:
                print(f"  {nm}: full timeseries at grid index (lat={li}, lon={oi})")

            for ti in range(nt):
                zval = float(da.isel(time=ti, lat=li, lon=oi).values)
                mm, dd, hh = _decode_time(da.time.values[ti])
                sval, existed = _source_value(root, fs, spec, mm, dd, hh,
                                              li, oi, sv)
                results["checked"] += 1
                stats["n"] += 1

                if not existed:
                    if np.isnan(zval):
                        results["nan_ok"] += 1
                    else:
                        results["mismatches"].append(
                            (nm, mm, dd, hh, li, oi,
                             f"source missing but zarr non-NaN value {zval}"))
                    continue

                src_nan = sval is None or np.isnan(sval)
                z_nan = np.isnan(zval)

                if src_nan and z_nan:
                    results["nan_value_ok"] += 1
                    stats["nan_both"] += 1
                    continue
                if z_nan and not src_nan:
                    results["mismatches"].append(
                        (nm, mm, dd, hh, li, oi,
                         f"zarr is NaN but source has real value {sval}"))
                    continue
                if src_nan and not z_nan:
                    results["mismatches"].append(
                        (nm, mm, dd, hh, li, oi,
                         f"source is NaN (fill) but zarr has value {zval}"))
                    continue

                stats["real"] += 1
                expected = _reverse_deaccum(zval, attrs)
                if np.isclose(expected, sval, atol=atol, rtol=rtol):
                    results["value_ok"] += 1
                else:
                    results["mismatches"].append(
                        (nm, mm, dd, hh, li, oi,
                         f"zarr={zval} (reversed={expected}) vs source={sval}"))

    return results


def check_deaccum_metadata(zarr_path, specs, var_groups, expect_no_deaccum,
                           verbose):
    """
    Independently re-derive the deaccumulation divisor from each
    variable's recorded source_variable/source_units and compare to
    what the store says was actually applied.

    This catches a class of bug the value round-trip check CANNOT:
    if reorganize_to_zarr.py computed the wrong divisor at conversion
    time, the written values and the recorded deaccumulation_seconds
    are wrong together, consistently - dividing back out by the same
    (wrong) recorded factor trivially reproduces the source value
    regardless of whether that factor was ever correct. Only an
    independent recomputation from the variable name / units can
    expose that.

    expect_no_deaccum : the store was deliberately built with
        --no-deaccumulate, so a variable that *could* have been
        converted but wasn't is expected, not a fault. Newer stores
        record this themselves in the 'deaccumulate_option' attr; this
        flag exists for stores written before that attr was added.
    """
    from reorganize_to_zarr import infer_deaccum_seconds
    issues, notes = [], []
    for nm in specs:
        if nm not in var_groups:
            continue
        da = _open_group(zarr_path, var_groups[nm])[nm]
        stored_flag = bool(da.attrs.get("deaccumulated", False))
        stored_secs = da.attrs.get("deaccumulation_seconds")
        sv = da.attrs.get("source_variable")
        src_units = da.attrs.get("source_units")
        recomputed = infer_deaccum_seconds(sv, src_units)

        # Did the store itself record that conversion was switched off?
        opt = da.attrs.get("deaccumulate_option")
        deliberately_off = (opt is False) or expect_no_deaccum

        if recomputed and not stored_flag:
            if deliberately_off:
                notes.append(
                    f"{nm}: stored RAW (accumulated) by choice. Values are "
                    f"'{src_units}' = J m-2 accumulated over "
                    f"{recomputed//3600}h, NOT W m-2. Divide by {recomputed} "
                    f"for average irradiance.")
            else:
                issues.append(
                    f"{nm}: stored deaccumulated=False but re-deriving from "
                    f"source_variable='{sv}' units='{src_units}' says it "
                    f"SHOULD be deaccumulated (divisor {recomputed}). If "
                    f"--no-deaccumulate was intended, pass "
                    f"--expect-no-deaccumulate.")
        elif stored_flag and not recomputed:
            issues.append(
                f"{nm}: stored deaccumulated=True but re-deriving from "
                f"'{sv}'/'{src_units}' says no conversion was warranted")
        elif recomputed and stored_secs != recomputed:
            issues.append(
                f"{nm}: stored deaccumulation_seconds={stored_secs} but "
                f"re-deriving from '{sv}'/'{src_units}' gives {recomputed}")
        elif verbose:
            if stored_flag:
                print(f"  {nm}: deaccumulated, divisor {stored_secs} - "
                      f"independently confirmed")
            else:
                print(f"  {nm}: not accumulated, no conversion expected - "
                      f"consistent")
    return issues, notes


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--zarr", required=True, help="path to the .zarr store")
    p.add_argument("--root", default=DEFAULT_ROOT,
                   help="original source root (local path or http(s) URL)")
    p.add_argument("--rad-var", default=None,
                   help="pass the same --rad-var used at conversion time, "
                        "if one was given")
    p.add_argument("--n-random", type=int, default=300,
                   help="random spot-check samples per variable")
    p.add_argument("--n-timeseries", type=int, default=5,
                   help="fixed grid points per variable to fully check "
                        "across every timestep")
    p.add_argument("--atol", type=float, default=1e-3)
    p.add_argument("--rtol", type=float, default=1e-5)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--expect-no-deaccumulate", action="store_true",
                   help="the store was built with --no-deaccumulate, so an "
                        "un-converted accumulated variable is expected "
                        "rather than a fault")
    p.add_argument("--vars", nargs="*", default=None,
                   choices=["2T", "SSRD", "WS10"])
    a = p.parse_args()

    specs = build_specs(a.rad_var)
    if a.vars:
        specs = {k: v for k, v in specs.items() if k in a.vars}
    fs = make_fs(a.root)

    print(f"Verifying {a.zarr}")
    print(f"against source: {a.root}")
    print(f"variables: {list(specs)}\n")

    print("=== Random spot checks ===")
    r1 = random_spot_checks(a.zarr, a.root, fs, specs, a.n_random,
                            a.atol, a.rtol, a.seed, verbose=True)
    print(f"  checked {r1['checked']}: {r1['value_ok']} real values matched, "
          f"{r1['nan_value_ok']} NaN-in-source matched, "
          f"{r1['nan_ok']} absent-file NaNs matched, "
          f"{len(r1['mismatches'])} mismatches")

    print("\n=== Full timeseries checks ===")
    r2 = full_timeseries_checks(a.zarr, a.root, fs, specs, a.n_timeseries,
                                a.atol, a.rtol, a.seed, verbose=True)
    print(f"  checked {r2['checked']}: {r2['value_ok']} real values matched, "
          f"{r2['nan_value_ok']} NaN-in-source matched, "
          f"{r2['nan_ok']} absent-file NaNs matched, "
          f"{len(r2['mismatches'])} mismatches")

    # NaN density is reported, not judged: a variable that is largely
    # NaN in BOTH store and source is faithfully converted, but it is
    # worth knowing about - it may be a legitimate masked domain, or it
    # may point at something upstream in the source archive.
    print("\n=== NaN density (store vs source agree; informational) ===")
    combined = {}
    for res in (r1, r2):
        for nm, s in res["per_var"].items():
            c = combined.setdefault(nm, {"n": 0, "nan_both": 0, "real": 0})
            for k in ("n", "nan_both", "real"):
                c[k] += s[k]
    for nm, s in combined.items():
        if s["n"]:
            pct = 100.0 * s["nan_both"] / s["n"]
            print(f"  {nm}: {s['nan_both']}/{s['n']} sampled points are NaN "
                  f"in both ({pct:.1f}%), {s['real']} real values compared")

    print("\n=== Deaccumulation metadata consistency ===")
    var_groups = discover_var_groups(a.zarr)
    meta_issues, meta_notes = check_deaccum_metadata(
        a.zarr, specs, var_groups, a.expect_no_deaccumulate, verbose=True)
    for note in meta_notes:
        print("  NOTE:", note)
    for issue in meta_issues:
        print("  ISSUE:", issue)

    all_mismatches = r1["mismatches"] + r2["mismatches"] + meta_issues
    if all_mismatches:
        print(f"\n*** {len(all_mismatches)} MISMATCHES FOUND ***")
        for m in all_mismatches[:30]:
            print(" ", m)
        if len(all_mismatches) > 30:
            print(f"  ... and {len(all_mismatches) - 30} more")
        sys.exit(1)
    else:
        total = r1["checked"] + r2["checked"]
        print(f"\nALL PASS - {total} total points checked, zero mismatches.")
        sys.exit(0)


if __name__ == "__main__":
    main()
