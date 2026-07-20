#!/usr/bin/env python3
"""Aggregate exp-D MM-GBSA/PBSA: 4 variants x 3 reps, Welch t, Spearman vs log(EC50)."""
import json, re, os, sys
import numpy as np
sys.path.insert(0, "/home/scroll/zzhang/semaglutide-in-silico/common/lib")
import stats as st
from scipy import stats as scst

BASE = os.path.dirname(os.path.abspath(__file__))
VARIANTS = ["no_linker", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]
EC50_PM = {"no_linker": 269.0, "gglu_1oeg": 4.8, "gglu_2oeg": 6.2, "gglu_3oeg": 27.7}  # Lau 2015 Table 3
REPS = [1, 2, 3]

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
                comp[name] = {"avg": float(m.group(2)), "std": float(m.group(3))}
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

res, summary = {}, {}
for v in VARIANTS:
    res[v] = {}
    for r in REPS:
        d = os.path.join(BASE, v, f"rep{r}")
        entry = {"final": parse_final(os.path.join(d, "FINAL_MMPBSA.dat"))}
        for m in ("GB", "PB"):
            x = delta_series(os.path.join(d, "FINAL_MMPBSA.csv"), m)
            q = st.summarize(x)
            entry[m] = {"mean": q["mean"], "std": q["std"], "n_eff": q["n_eff"],
                        "sem": q["std"] / np.sqrt(q["n_eff"])}
        res[v][f"rep{r}"] = entry
    summary[v] = {}
    for m in ("GB", "PB"):
        means = np.array([res[v][f"rep{r}"][m]["mean"] for r in REPS])
        summary[v][m] = {"replica_means": means.tolist(), "mean": float(means.mean()),
                         "sd": float(means.std(ddof=1)), "sem": float(means.std(ddof=1) / np.sqrt(3))}

json.dump({"per_rep": res, "summary": summary}, open(os.path.join(BASE, "aggregated.json"), "w"), indent=1)

print(f"{'variant':12s} {'dG_GB (3rep)':>18s} {'dG_PB (3rep)':>18s}   EC50(pM)")
for v in VARIANTS:
    print(f"{v:12s} {summary[v]['GB']['mean']:9.2f}±{summary[v]['GB']['sd']:5.2f}   "
          f"{summary[v]['PB']['mean']:9.2f}±{summary[v]['PB']['sd']:5.2f}   {EC50_PM[v]}")

# pairwise Welch t vs no_linker
print("\npairwise Welch t (variant - no_linker):")
for v in VARIANTS[1:]:
    for m in ("GB", "PB"):
        a = summary[v][m]["replica_means"]; b = summary["no_linker"][m]["replica_means"]
        t, p = scst.ttest_ind(a, b, equal_var=False)
        dd = np.mean(a) - np.mean(b)
        print(f"  {v:12s} {m}: ddG={dd:+7.2f}  t={t:+.2f}  p={p:.4f}")

# Spearman: variant means (n=4) and all replica points (n=12)
for m in ("GB", "PB"):
    d_mean = [summary[v][m]["mean"] for v in VARIANTS]
    e_mean = [np.log10(EC50_PM[v]) for v in VARIANTS]
    R4, p4 = scst.spearmanr(d_mean, e_mean)
    d_all = [res[v][f"rep{r}"][m]["mean"] for v in VARIANTS for r in REPS]
    e_all = [np.log10(EC50_PM[v]) for v in VARIANTS for r in REPS]
    R12, p12 = scst.spearmanr(d_all, e_all)
    print(f"\nSpearman dG_{m} vs log10(EC50): variant-mean R={R4:+.3f} (n=4, p={p4:.3f}) | all-replica R={R12:+.3f} (n=12, p={p12:.4f})")
    summary.setdefault("spearman", {})[m] = {"R_n4": float(R4), "p_n4": float(p4), "R_n12": float(R12), "p_n12": float(p12)}

json.dump(summary, open(os.path.join(BASE, "summary.json"), "w"), indent=1)
