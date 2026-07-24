#!/usr/bin/env python3
"""Aggregate chain-series MM-GBSA/PBSA: 5 diacids x 3 reps, pairwise interactions, decomp.
True total carbon numbers: c12=12, c14=14, c16=16, c20=20, c22=22.
Note: legacy "c18_diacid" (exp-C wave-0) is physically the SAME C20 diacid molecule;
its results (GB -89.69+-4.02, PB -55.64+-7.35) serve as cross-validation for c20_diacid.
"""
import json, re, os, sys
import numpy as np
sys.path.insert(0, "/home/scroll/zzhang/semaglutide-in-silico/common/lib")
import stats as st

BASE = os.path.dirname(os.path.abspath(__file__))
SYSTEMS = ["c12_diacid", "c14_diacid", "c16_diacid", "c20_diacid", "c22_diacid"]
CARBONS = {"c12_diacid": 12, "c14_diacid": 14, "c16_diacid": 16, "c20_diacid": 20, "c22_diacid": 22}
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

def parse_pair(path, tag):
    hdr = open(path).readline().split()
    data = np.loadtxt(path)
    gv = lambda k: [i for i, x in enumerate(hdr) if x == k][0]
    vdw = data[:, gv(f"{tag}_ab[vdw]")] - data[:, gv(f"{tag}_a[vdw]")] - data[:, gv(f"{tag}_b[vdw]")]
    ele = data[:, gv(f"{tag}_ab[elec]")] - data[:, gv(f"{tag}_a[elec]")] - data[:, gv(f"{tag}_b[elec]")]
    return {"vdw": float(vdw.mean()), "ele": float(ele.mean())}

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
                sections[method][f[0].strip()] = {"vdw": float(f[5]), "ele": float(f[8]),
                    "polar": float(f[11]), "nonpolar": float(f[14]), "total": float(f[17])}
            except (ValueError, IndexError):
                continue
    return sections

res, summary = {}, {}
for s in SYSTEMS:
    res[s] = {}
    for r in REPS:
        d = os.path.join(BASE, s, f"rep{r}")
        entry = {"final": parse_final(os.path.join(d, "FINAL_MMPBSA.dat")),
                 "decomp": parse_decomp(os.path.join(d, "FINAL_DECOMP.dat")),
                 "pw_distal": parse_pair(os.path.join(d, "pw_distal_anchor.dat"), "e"),
                 "pw_prox": parse_pair(os.path.join(d, "pw_prox_408.dat"), "p")}
        for m in ("GB", "PB"):
            x = delta_series(os.path.join(d, "FINAL_MMPBSA.csv"), m)
            q = st.summarize(x)
            entry[m] = {"mean": q["mean"], "std": q["std"], "n_eff": q["n_eff"],
                        "sem": q["std"] / np.sqrt(q["n_eff"])}
        res[s][f"rep{r}"] = entry
    summary[s] = {"carbons": CARBONS[s]}
    for m in ("GB", "PB"):
        means = np.array([res[s][f"rep{r}"][m]["mean"] for r in REPS])
        summary[s][m] = {"replica_means": means.tolist(), "mean": float(means.mean()),
                         "sd": float(means.std(ddof=1)), "sem": float(means.std(ddof=1) / np.sqrt(3))}
    for pw in ("pw_distal", "pw_prox"):
        ele = np.array([res[s][f"rep{r}"][pw]["ele"] for r in REPS])
        vdw = np.array([res[s][f"rep{r}"][pw]["vdw"] for r in REPS])
        summary[s][pw] = {"ele_mean": float(ele.mean()), "ele_sd": float(ele.std(ddof=1)),
                          "vdw_mean": float(vdw.mean()), "vdw_sd": float(vdw.std(ddof=1))}

json.dump({"per_rep": res, "summary": summary}, open(os.path.join(BASE, "aggregated.json"), "w"), indent=1)

print(f"{'system':12s} {'C':>3s} {'dG_GB mean±SD':>16s} {'dG_PB mean±SD':>16s}  replicas(GB)")
for s in SYSTEMS:
    q = summary[s]
    print(f"{s:12s} {q['carbons']:3d} {q['GB']['mean']:9.2f}±{q['GB']['sd']:5.2f}   {q['PB']['mean']:9.2f}±{q['PB']['sd']:5.2f}   "
          + " ".join(f"{x:.1f}" for x in q['GB']['replica_means']))

print("\npairwise (bare MM ele, 3-rep mean±SD):")
for s in SYSTEMS:
    q = summary[s]
    print(f"{s:12s} distalCOO-ARG346/483 ele {q['pw_distal']['ele_mean']:8.2f}±{q['pw_distal']['ele_sd']:5.2f} | "
          f"proxCOO-ARG408 ele {q['pw_prox']['ele_mean']:8.2f}±{q['pw_prox']['ele_sd']:6.2f}")

print("\ndechomp key residues (TDC total, 3-rep mean; GB):")
keys = ["ARG 346", "ARG 408", "ARG 483", "FAH 583"]
print(f"{'system':12s} " + " ".join(f"{k:>12s}" for k in keys))
for s in SYSTEMS:
    line = f"{s:12s} "
    for k in keys:
        vals = [res[s][f"rep{r}"]["decomp"].get("GB", {}).get(k, {}).get("total") for r in REPS]
        vals = [v for v in vals if v is not None]
        line += f"{np.mean(vals):12.2f} " if vals else "       N/A  "
    print(line)

# cross-check c20_diacid vs legacy c18_diacid
print("\ncross-check: c20_diacid GB {0:.2f}±{1:.2f} / PB {2:.2f}±{3:.2f}  vs  legacy c18_diacid GB -89.69±4.02 / PB -55.64±7.35".format(
    summary["c20_diacid"]["GB"]["mean"], summary["c20_diacid"]["GB"]["sd"],
    summary["c20_diacid"]["PB"]["mean"], summary["c20_diacid"]["PB"]["sd"]))
