#!/usr/bin/env python3
"""Aggregate exp-C MM-GBSA/PBSA results: per-replica stats, 3-replica means, ddG, decomp."""
import json, re, os, sys
import numpy as np
sys.path.insert(0, "/home/scroll/zzhang/semaglutide-in-silico/common/lib")
import stats as st

BASE = os.path.dirname(os.path.abspath(__file__))
SYSTEMS = ["c18_monoacid", "c18_diacid"]
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

def parse_decomp(path):
    lines = open(path).read().splitlines()
    sections, method, in_tdc = {}, None, False
    for ln in lines:
        if "Energy Decomposition Analysis" in ln:
            method = "GB" if "Generalized Born" in ln else ("PB" if "Poisson" in ln else method)
            sections.setdefault(method, {}); in_tdc = False; continue
        if ln.startswith("Total Energy Decomposition:"):
            in_tdc = True; continue
        if ln.startswith(("Sidechain Energy Decomposition:", "Backbone Energy Decomposition:")):
            in_tdc = False; continue
        if in_tdc and re.match(r"^[A-Z0-9]{3} \d+,", ln):
            f = ln.split(",")
            try:
                sections[method][f[0].strip()] = {
                    "vdw": float(f[5]), "ele": float(f[8]), "polar": float(f[11]),
                    "nonpolar": float(f[14]), "total": float(f[17]), "total_std": float(f[18])}
            except (ValueError, IndexError):
                continue
    return sections

res = {}
for s in SYSTEMS:
    res[s] = {}
    for r in REPS:
        d = os.path.join(BASE, s, f"rep{r}")
        entry = {"final": parse_final(os.path.join(d, "FINAL_MMPBSA.dat")),
                 "decomp": parse_decomp(os.path.join(d, "FINAL_DECOMP.dat"))}
        for m in ("GB", "PB"):
            x = delta_series(os.path.join(d, "FINAL_MMPBSA.csv"), m)
            q = st.summarize(x)
            entry[m] = {"mean": q["mean"], "std": q["std"], "n_eff": q["n_eff"],
                        "sem": q["std"] / np.sqrt(q["n_eff"])}
        res[s][f"rep{r}"] = entry

json.dump(res, open(os.path.join(BASE, "aggregated.json"), "w"), indent=1)

# ---- replica-level summary ----
print(f"{'system':14s} {'rep':4s} {'dG_GB':>10s} {'SEM':>6s} {'dG_PB':>10s} {'SEM':>6s}")
summary = {}
for s in SYSTEMS:
    for r in REPS:
        e = res[s][f"rep{r}"]
        print(f"{s:14s} rep{r}  {e['GB']['mean']:10.2f} {e['GB']['sem']:6.2f} {e['PB']['mean']:10.2f} {e['PB']['sem']:6.2f}")
    for m in ("GB", "PB"):
        means = np.array([res[s][f"rep{r}"][m]["mean"] for r in REPS])
        summary.setdefault(s, {})[m] = {"replica_means": means.tolist(),
                                        "mean": float(means.mean()), "sd": float(means.std(ddof=1)),
                                        "sem": float(means.std(ddof=1) / np.sqrt(len(means)))}
print()
for m in ("GB", "PB"):
    a = summary["c18_diacid"][m]; b = summary["c18_monoacid"][m]
    dd = a["mean"] - b["mean"]
    dd_sem = np.sqrt(a["sem"] ** 2 + b["sem"] ** 2)
    # Welch t on replica means (n=3 vs 3)
    from scipy import stats as scst
    t, p = scst.ttest_ind(a["replica_means"], b["replica_means"], equal_var=False)
    print(f"ddG({m}) = diacid - monoacid = {dd:+.2f} +- {dd_sem:.2f} kcal/mol (replica-mean Welch t={t:.2f}, p={p:.4f})")
    summary.setdefault("ddG", {})[m] = {"ddG": float(dd), "sem": float(dd_sem), "t": float(t), "p": float(p)}

json.dump(summary, open(os.path.join(BASE, "summary.json"), "w"), indent=1)

# ---- decomp comparison (3-replica means) ----
print("\nkey residue decomp (TDC total, 3-replica mean; GB / PB):")
keys = ["SER 340", "ARG 346", "LYS 349", "ARG 408", "TYR 409", "ARG 482", "ARG 483", "FAH 583"]
for k in keys:
    line = f"  {k:9s}"
    for s in SYSTEMS:
        for m in ("GB", "PB"):
            vals = [res[s][f"rep{r}"]["decomp"].get(m, {}).get(k, {}).get("total") for r in REPS]
            vals = [v for v in vals if v is not None]
            if vals:
                line += f" | {s.split('_')[1][:4]} {m} {np.mean(vals):7.2f}"
    print(line)
