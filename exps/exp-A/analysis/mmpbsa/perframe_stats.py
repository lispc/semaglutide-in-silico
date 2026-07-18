#!/usr/bin/env python3
"""Per-frame statistics for exp-A MM-PBSA: n_eff-aware errors, paired t-test, clash correction."""
import sys, numpy as np
sys.path.insert(0, "/home/scroll/zzhang/semaglutide-in-silico/common/lib")
import stats

BASE = "/home/scroll/zzhang/semaglutide-in-silico/exps/exp-A/analysis/mmpbsa"

def delta_series(csv_path, method):
    """Return per-frame DELTA TOTAL array from FINAL_MMPBSA.csv for method GB/PB."""
    lines = open(csv_path).read().splitlines()
    start = None
    for i, ln in enumerate(lines):
        if method == "GB" and ln.startswith("GENERALIZED BORN"):
            start = i
        if method == "PB" and ln.startswith("POISSON BOLTZMANN"):
            start = i
    # find the Deltas block header after `start`
    for i in range(start, len(lines)):
        if lines[i].startswith("Frame #") and "DELTA TOTAL" in lines[i]:
            hdr = i
            break
    ncol = lines[hdr].split(",").index("DELTA TOTAL")
    vals = []
    for ln in lines[hdr + 1:]:
        f = ln.split(",")
        if len(f) <= ncol or not f[0].strip().isdigit():
            break
        vals.append(float(f[ncol]))
    return np.array(vals)

def pair_series(path):
    vdws, eles = [], []
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.split()
        vdws.append(float(f[6])); eles.append(float(f[7]))
    return np.array(vdws) + np.array(eles)

out = {}
for m in ("GB", "PB"):
    s = {}
    for sys_ in ("wt", "aib8"):
        x = delta_series(f"{BASE}/{sys_}/FINAL_MMPBSA.csv", m)
        s[sys_] = stats.summarize(x, f"{sys_}_{m}")
        s[sys_]["series"] = x
    d = s["aib8"]["series"] - s["wt"]["series"]
    t = stats.correlated_t_test(s["aib8"]["series"], s["wt"]["series"])
    dsum = stats.summarize(d, f"ddG_{m}")
    out[m] = {"wt": {k: v for k, v in s["wt"].items() if k != "series"},
              "aib8": {k: v for k, v in s["aib8"].items() if k != "series"},
              "ddG_mean": float(d.mean()), "ddG_std": float(d.std(ddof=1)),
              "ddG_n_eff": dsum["n_eff"], "ddG_sem": float(d.std(ddof=1) / np.sqrt(dsum["n_eff"])),
              "ttest": t}
    print(f"== {m} ==")
    for sys_ in ("wt", "aib8"):
        q = s[sys_]
        print(f"  {sys_:5s} mean={q['mean']:8.2f} std={q['std']:6.2f} median={q['median']:8.2f} "
              f"IQR=[{q['q25']:.1f},{q['q75']:.1f}] tau={q['tau_frames']:.1f} n_eff={q['n_eff']:.1f}")
    print(f"  ddG = {d.mean():+.2f} +- {d.std(ddof=1):.2f} (n={len(d)}), "
          f"SEM(n_eff)={out[m]['ddG_sem']:.2f}, paired t={t['t']:.2f}, p={t['p']:.3f}, n_eff={t['n_eff']:.1f}")

# clash-corrected ddG (subtract pair MM interaction difference; also correct per-frame for GB)
pw = pair_series(f"{BASE}/wt/pair_658_750.dat")
pa = pair_series(f"{BASE}/aib8/pair_658_750.dat")
dd_clash = (pa - pw)
print(f"\nclash pair ddE per-frame: mean={dd_clash.mean():+.2f} std={dd_clash.std(ddof=1):.2f}")
for m in ("GB", "PB"):
    print(f"  clash-corrected ddG({m}) = {out[m]['ddG_mean'] - dd_clash.mean():+.2f} kcal/mol")

import json
json.dump({m: {k: v for k, v in out[m].items()} for m in out},
          open(f"{BASE}/perframe_stats.json", "w"), indent=1, default=float)
