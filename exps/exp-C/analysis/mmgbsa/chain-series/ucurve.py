#!/usr/bin/env python3
"""Chain-series 9-point summary + U-curve plot (GB & PB panels, x = true total carbons)."""
import json, re, os, sys
import numpy as np
sys.path.insert(0, "/home/scroll/zzhang/semaglutide-in-silico/common/lib")
import stats as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
MMGBSA = os.path.dirname(BASE)  # analysis/mmgbsa (legacy c18 systems live here)

def parse_final(path):
    txt = open(path).read()
    out = {}
    for method, key in [("GENERALIZED BORN:", "GB"), ("POISSON BOLTZMANN:", "PB")]:
        i = txt.find(method)
        if i < 0:
            continue
        j = txt.find("Differences (Complex - Receptor - Ligand):", i)
        comp = {}
        for line in txt[j:j + 2200].splitlines():
            m = re.match(r"([A-Z][A-Za-z0-9 ]+?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
            if m:
                name = m.group(1).strip()
                if name.startswith(("Complex", "Receptor", "Ligand", "Energy Component", "---")):
                    continue
                comp[name] = {"avg": float(m.group(2))}
                if name == "DELTA TOTAL":
                    break
        out[key] = comp
    return out

def delta_series(csv_path, method):
    lines = open(csv_path).read().splitlines()
    start = 0
    for i, ln in enumerate(lines):
        if (method == "GB" and ln.startswith("GENERALIZED BORN")) or \
           (method == "PB" and ln.startswith("POISSON BOLTZMANN")):
            start = i
    hdr = None
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

def system_stats(sysdir, reps):
    out = {}
    for m in ("GB", "PB"):
        means = []
        for r in reps:
            x = delta_series(os.path.join(sysdir, f"rep{r}", "FINAL_MMPBSA.csv"), m)
            means.append(float(x.mean()))
        means = np.array(means)
        out[m] = {"replica_means": means.tolist(), "mean": float(means.mean()),
                  "sd": float(means.std(ddof=1)) if len(means) > 1 else 0.0, "n": len(means)}
    return out

# --- collect all points ---
points = []  # (series, carbon, system_label, stats)
DI = [("c12_diacid", 12, [1,2,3]), ("c14_diacid", 14, [1,2,3]), ("c16_diacid", 16, [1,2,3]),
      ("c18true_diacid", 18, [1,2]), ("c20_diacid", 20, [1,2,3]), ("c22_diacid", 22, [1,2,3])]
for s, c, reps in DI:
    points.append(("diacid", c, s, system_stats(os.path.join(BASE, s), reps)))
# legacy c18_diacid = C20 (physically same molecule as c20_diacid)
points.append(("diacid", 20, "c18_diacid(legacy=C20)", system_stats(os.path.join(MMGBSA, "c18_diacid"), [1,2,3])))
# monoacids
points.append(("monoacid", 16, "c16_monoacid", system_stats(os.path.join(BASE, "c16_monoacid"), [1,2,3])))
points.append(("monoacid", 19, "c18_monoacid(legacy=C19)", system_stats(os.path.join(MMGBSA, "c18_monoacid"), [1,2,3])))

json.dump({p[2]: {"series": p[0], "carbons": p[1], **p[3]} for p in points},
          open(os.path.join(BASE, "ucurve_points.json"), "w"), indent=1)

print(f"{'series':9s} {'C':>3s} {'system':26s} {'n':>2s} {'dG_GB':>15s} {'dG_PB':>15s}")
for series, c, label, q in sorted(points, key=lambda p: (p[0], p[1])):
    print(f"{series:9s} {c:3d} {label:26s} {q['GB']['n']:2d} "
          f"{q['GB']['mean']:8.2f}±{q['GB']['sd']:5.2f} {q['PB']['mean']:8.2f}±{q['PB']['sd']:5.2f}")

# --- plot ---
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=False)
colors = {"diacid": "#1f77b4", "monoacid": "#d62728"}
for ax, m, title in [(axes[0], "GB", "MM-GBSA (igb=5)"), (axes[1], "PB", "MM-PBSA (inp=2)")]:
    for series, marker in [("diacid", "o"), ("monoacid", "s")]:
        pts = sorted([p for p in points if p[0] == series], key=lambda p: p[1])
        # legacy C20 duplicates offset for visibility
        xs = [p[1] + (-0.15 if "legacy" in p[2] else (0.15 if p[2] == "c20_diacid" else 0)) for p in pts]
        ys = [p[3][m]["mean"] for p in pts]
        es = [p[3][m]["sd"] for p in pts]
        ax.errorbar(xs, ys, yerr=es, fmt=marker + "-", color=colors[series], capsize=3,
                    ms=6, lw=1.2, label={"diacid": "diacid series", "monoacid": "monoacid"}[series])
        for x, y, p in zip(xs, ys, pts):
            note = "*" if p[3][m]["n"] < 3 else ""
            ax.annotate(f"C{p[1]}{note}", (x, y), textcoords="offset points", xytext=(0, 7),
                        ha="center", fontsize=8, color=colors[series])
    ax.set_xlabel("true total carbon number")
    ax.set_ylabel("ΔG_bind (kcal/mol)")
    ax.set_title(f"{title} — FA3 chain-length series")
    ax.set_xticks([12, 14, 16, 18, 20, 22])
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="lower left" if m == "GB" else "upper left")
fig.suptitle("exp-C chain series @ HSA FA3 (last 50 ns ×3 reps; *=n=2, rep3 pending)", fontsize=9)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(os.path.join(BASE, "U_curve.png"), dpi=160)
print("\nsaved U_curve.png")
