#!/usr/bin/env python3
"""Parse exp-A MM-PBSA/GBSA outputs (FINAL_MMPBSA.dat + FINAL_DECOMP.dat) for WT and Aib8.

Outputs a JSON summary used to fill RESULTS.md:
- GB / PB delta components (mean, std) for both systems
- per-residue TDC (GB and PB) for printed residues
- clash pair (Lys696<->Phe28) correction from pair_658_750.dat
"""
import json, re, os, sys

BASE = os.path.dirname(os.path.abspath(__file__))
SYSTEMS = ["wt", "aib8"]

def parse_final(path):
    """Extract GB and PB 'Differences' blocks from FINAL_MMPBSA.dat."""
    txt = open(path).read()
    out = {}
    for method, key in [("GENERALIZED BORN:", "GB"), ("POISSON BOLTZMANN:", "PB")]:
        i = txt.find(method)
        if i < 0:
            continue
        j = txt.find("Differences (Complex - Receptor - Ligand):", i)
        block = txt[j:j + 2200]
        comp = {}
        for line in block.splitlines():
            m = re.match(r"([A-Z][A-Za-z0-9 ]+?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", line)
            if m:
                name = m.group(1).strip()
                if name.startswith(("Complex", "Receptor", "Ligand", "Energy Component", "---")):
                    continue
                comp[name] = {"avg": float(m.group(2)), "std": float(m.group(3))}
                if name == "DELTA TOTAL":
                    break  # stop before the next section's blocks overwrite keys
        out[key] = comp
    return out

def parse_decomp(path):
    """Parse FINAL_DECOMP.dat: per-residue rows for GB and PB sections (TDC block)."""
    lines = open(path).read().splitlines()
    sections = {}  # 'GB'/'PB' -> {residue_label: {term: avg}}
    method = None
    in_tdc = False
    for ln in lines:
        if "Energy Decomposition Analysis" in ln:
            method = "GB" if "Generalized Born" in ln else ("PB" if "Poisson" in ln else method)
            sections.setdefault(method, {})
            in_tdc = False
            continue
        if ln.startswith("Total Energy Decomposition:"):
            in_tdc = True
            continue
        if ln.startswith(("Sidechain Energy Decomposition:", "Backbone Energy Decomposition:")):
            in_tdc = False
            continue
        if in_tdc and re.match(r"^[A-Z0-9]{3} \d+,", ln):
            f = ln.split(",")
            # columns: name, location, Internal(3), vdW(3), Ele(3), Polar(3), NonPolar(3), TOTAL(3)
            try:
                vals = {
                    "vdw": float(f[5]), "ele": float(f[8]),
                    "polar": float(f[11]), "nonpolar": float(f[14]),
                    "total": float(f[17]), "total_std": float(f[18]),
                }
            except (ValueError, IndexError):
                continue
            sections[method][f[0].strip() + "|" + f[1].strip()] = vals
    return sections

def parse_pair(path):
    """Mean pairwise vdW/elec for the Lys696-Phe28 clash pair."""
    vdws, eles = [], []
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.split()
        vdws.append(float(f[6])); eles.append(float(f[7]))
    n = len(vdws)
    m = sum(vdws) / n
    sd = (sum(x * x for x in vdws) / n - m * m) ** 0.5
    return {"vdw_avg": m, "vdw_std": sd, "ele_avg": sum(eles) / n, "n": n}

res = {}
for s in SYSTEMS:
    d = os.path.join(BASE, s)
    res[s] = {
        "final": parse_final(os.path.join(d, "FINAL_MMPBSA.dat")),
        "decomp": parse_decomp(os.path.join(d, "FINAL_DECOMP.dat")),
        "pair": parse_pair(os.path.join(d, "pair_658_750.dat")),
    }

json.dump(res, open(os.path.join(BASE, "parsed_results.json"), "w"), indent=1)

# human-readable quick summary
for s in SYSTEMS:
    for m in ("GB", "PB"):
        dt = res[s]["final"].get(m, {}).get("DELTA TOTAL", {})
        print(f"{s:5s} {m} DELTA TOTAL: {dt.get('avg'):>9.2f} +- {dt.get('std'):.2f}")
for m in ("GB", "PB"):
    a = res["aib8"]["final"][m]["DELTA TOTAL"]["avg"]
    w = res["wt"]["final"][m]["DELTA TOTAL"]["avg"]
    sa = res["aib8"]["final"][m]["DELTA TOTAL"]["std"]
    sw = res["wt"]["final"][m]["DELTA TOTAL"]["std"]
    print(f"ddG({m}) = {a - w:+.2f} +- {(sa**2 + sw**2)**0.5:.2f} kcal/mol (std-propagated)")
pw, pa = res["wt"]["pair"], res["aib8"]["pair"]
clash_dd = (pa["vdw_avg"] + pa["ele_avg"]) - (pw["vdw_avg"] + pw["ele_avg"])
print(f"clash pair ddE(MM) = {clash_dd:+.2f} kcal/mol (Aib8 - WT)")
