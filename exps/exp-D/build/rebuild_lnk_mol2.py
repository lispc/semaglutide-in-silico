#!/usr/bin/env python3
"""Rebuild LNK mol2 inputs for exp-D restart (2026-07-17).

The original tleap/lnk_*_pos.mol2 inputs are structurally broken:
- 13-17 orphan H atoms per file (leftover from a name-based strip whose
  hardcoded DELETE list only matched the no_linker/gglu atom layout)
- gglu_1oeg: mid-chain atoms C32, O33, C35 deleted -> chain in fragments
- gglu_1oeg/2oeg/3oeg: the NME cap fragment (last 4 heavy atoms) was never
  deleted and floats ~7 A from the linker
- no BOND table; all internal amides mistyped c3/o (free rotation)

This script redoes the strip structure-aware from build/{v}_bcc.mol2
(same RDKit conformers the original pos mol2 came from, so kept-atom
geometry is identical) + bond table from build/{v}_raw.mol2:
- delete ACE (C1,C2,O3), Lys backbone (N4,C5), Lys sidechain (C6-C9),
  the sidechain NZ (N10; the protein's Lys26 NZ takes its place), and the
  backbone NME cap (last 4 heavy atoms, per-variant names)
- delete all H bonded to deleted heavy atoms
- assign proper GAFF2 types from the RDKit bond orders (c/o for carbonyls,
  os for ether O, n for amide N, n3 for the gGlu amine, c3/hc/hn)
- position via the original transform (N10->C11 vector aligned to the
  peptide CE->NZ direction, N10 placed at NZ position)
- charges kept verbatim from the bcc file (same charge model the original
  production used; per-variant net-charge offset is neutralized by tleap)

Output: tleap/lnk_{variant}_clean.mol2 (heavy atoms first, then H)
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TLEAP = os.path.normpath(os.path.join(HERE, "..", "tleap"))

# Lys26 NZ / CE positions in ecd_pep.pdb (from position_lnk.py)
NZ_POS = np.array([-7.571, -0.448, 1.906])
CE_POS = np.array([-6.269, -0.196, 1.405])
CE_TO_NZ = NZ_POS - CE_POS
CE_TO_NZ_DIR = CE_TO_NZ / np.linalg.norm(CE_TO_NZ)

DELETE_COMMON = ["C1", "C2", "O3", "N4", "C5", "C6", "C7", "C8", "C9", "N10"]
VARIANTS = {
    # name: NME cap heavy atoms (backbone C(=O)-NH-CH3, last 4 heavy atoms)
    "no_linker": ["C32", "O33", "N34", "C35"],
    "gglu_1oeg": ["C48", "O49", "N50", "C51"],
    "gglu_2oeg": ["C58", "O59", "N60", "C61"],
    "gglu_3oeg": ["C68", "O69", "N70", "C71"],
}

def parse_mol2(path):
    atoms, bonds, sec = {}, [], None
    for line in open(path):
        s = line.strip()
        if s.startswith("@<TRIPOS>"):
            sec = s; continue
        if sec == "@<TRIPOS>ATOM" and s:
            p = s.split()
            atoms[int(p[0])] = dict(name=p[1],
                                    xyz=np.array([float(p[2]), float(p[3]), float(p[4])]),
                                    type=p[5], charge=p[8])
        elif sec == "@<TRIPOS>BOND" and s:
            p = s.split()
            bonds.append((int(p[1]), int(p[2]), p[3]))
    return atoms, bonds

def gaff2_type(name, heavy_bonds, bond_orders):
    """Assign GAFF2 type from local chemistry. heavy_bonds: dict name->set of
    bonded heavy names; bond_orders: dict frozenset({a,b}) -> order ('1'/'2')."""
    elem = name[0]
    if elem == "H":
        return "hc"  # H on N reassigned by caller context
    if elem == "C":
        for nb in heavy_bonds[name]:
            if nb.startswith("O") and bond_orders.get(frozenset([name, nb])) == "2":
                return "c"  # carbonyl / carboxyl C
        return "c3"
    if elem == "O":
        for nb in heavy_bonds[name]:
            if bond_orders.get(frozenset([name, nb])) == "2":
                return "o"  # double-bonded O
        # single-bonded O: carboxylate (neighbor is 'c') or ether
        for nb in heavy_bonds[name]:
            if any(bond_orders.get(frozenset([nb, x])) == "2" for x in heavy_bonds[nb] if x.startswith("O")):
                return "o"  # carboxylate O(-)
        return "os"
    if elem == "N":
        for nb in heavy_bonds[name]:
            if any(bond_orders.get(frozenset([nb, x])) == "2" and x.startswith("O") for x in heavy_bonds[nb]):
                return "n"  # amide N
        return "n3"  # amine N (gGlu alpha-NH2)
    return "du"

def rebuild(v):
    bcc_atoms, _ = parse_mol2(os.path.join(HERE, f"{v}_bcc.mol2"))
    raw_atoms, raw_bonds = parse_mol2(os.path.join(HERE, f"{v}_raw.mol2"))
    assert set(bcc_atoms) == set(raw_atoms), f"{v}: bcc/raw atom id mismatch"
    for i in bcc_atoms:
        assert bcc_atoms[i]["name"] == raw_atoms[i]["name"], f"{v}: name mismatch at {i}"

    heavy = {i: a for i, a in bcc_atoms.items() if not a["name"].startswith("H")}
    hs = {i: a for i, a in bcc_atoms.items() if a["name"].startswith("H")}

    # Heavy-atom bond graph + bond orders (from raw mol2, written by RDKit)
    heavy_bonds = {a["name"]: set() for a in heavy.values()}
    bond_orders = {}
    h_to_heavy = {}
    idname = {i: a["name"] for i, a in bcc_atoms.items()}
    for a1, a2, bt in raw_bonds:
        n1, n2 = idname[a1], idname[a2]
        if a1 in heavy and a2 in heavy:
            heavy_bonds[n1].add(n2); heavy_bonds[n2].add(n1)
            bond_orders[frozenset([n1, n2])] = bt
        elif a1 in hs and a2 in heavy:
            h_to_heavy[n1] = n2
        elif a2 in hs and a1 in heavy:
            h_to_heavy[n2] = n1

    # Verify NME cap pattern: C(=O)-N-C at the chain backbone end
    nme = VARIANTS[v]
    assert [n[0] for n in nme] == ["C", "O", "N", "C"], f"{v}: NME pattern {nme}"
    assert bond_orders.get(frozenset([nme[0], nme[1]])) == "2", f"{v}: {nme[0]}={nme[1]} not double"

    delete_names = set(DELETE_COMMON + nme)
    delete_ids = {i for i, a in bcc_atoms.items() if a["name"] in delete_names}
    delete_ids |= {i for i, a in hs.items() if h_to_heavy.get(a["name"]) in delete_names}

    keep = [(i, a) for i, a in bcc_atoms.items() if i not in delete_ids]
    keep_heavy = [(i, a) for i, a in keep if not a["name"].startswith("H")]
    keep_h = [(i, a) for i, a in keep if a["name"].startswith("H")]

    # Sanity: kept heavy graph is a single connected component
    keep_names = {a["name"] for _, a in keep_heavy}
    adj = {n: heavy_bonds[n] & keep_names for n in keep_names}
    seen, stack = set(), [next(iter(keep_names))]
    while stack:
        u = stack.pop()
        if u in seen: continue
        seen.add(u); stack.extend(adj[u] - seen)
    assert seen == keep_names, f"{v}: kept graph not connected ({len(seen)}/{len(keep_names)})"

    # Types (H on N -> hn, H on C -> hc)
    types = {}
    for _, a in keep_heavy:
        types[a["name"]] = gaff2_type(a["name"], heavy_bonds, bond_orders)
    for _, a in keep_h:
        parent = h_to_heavy[a["name"]]
        types[a["name"]] = "hn" if parent.startswith("N") else "hc"

    # Position: original transform (N10 -> NZ_POS, N10->C11 along CE->NZ)
    n_pos = bcc_atoms[[i for i, a in bcc_atoms.items() if a["name"] == "N10"][0]]["xyz"]
    c_pos = bcc_atoms[[i for i, a in bcc_atoms.items() if a["name"] == "C11"][0]]["xyz"]
    v1 = (c_pos - n_pos) / np.linalg.norm(c_pos - n_pos)
    v2 = CE_TO_NZ_DIR
    cos_a = float(np.dot(v1, v2))
    if abs(cos_a - 1) < 1e-8:
        R = np.eye(3)
    elif abs(cos_a + 1) < 1e-8:
        R = -np.eye(3)
    else:
        k = np.cross(v1, v2); k /= np.linalg.norm(k)
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        R = np.eye(3) + np.sin(np.arccos(cos_a)) * K + (1 - cos_a) * K @ K

    body = keep_heavy + keep_h  # heavy first, then H
    out_atoms, old2new = [], {}
    for new_i, (old_i, a) in enumerate(body, 1):
        old2new[old_i] = new_i
        out_atoms.append(dict(name=a["name"], xyz=R @ (a["xyz"] - n_pos) + NZ_POS,
                              type=types[a["name"]], charge=a["charge"]))
    out_bonds = [(old2new[a1], old2new[a2], bt) for a1, a2, bt in raw_bonds
                 if a1 in old2new and a2 in old2new]

    # Geometry sanity: all bonds short; NZ->C11 distance preserved
    for a1, a2, _ in out_bonds:
        d = np.linalg.norm(out_atoms[a1 - 1]["xyz"] - out_atoms[a2 - 1]["xyz"])
        assert d < 1.7, f"{v}: long bond {out_atoms[a1-1]['name']}-{out_atoms[a2-1]['name']} {d:.2f}"
    nz_c = np.linalg.norm(out_atoms[0]["xyz"] - NZ_POS)  # atom order: C11 first?
    netq = sum(float(a["charge"]) for a in out_atoms)

    out = os.path.join(TLEAP, f"lnk_{v}_clean.mol2")
    with open(out, "w") as f:
        f.write("@<TRIPOS>MOLECULE\nLNK\n")
        f.write(f" {len(out_atoms):5d} {len(out_bonds):5d}     0     0     0\nSMALL\nGAFF2\n")
        f.write(f"rebuilt {v} linker-C18 (structure-aware strip, GAFF2 types, bonds)\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for i, a in enumerate(out_atoms, 1):
            f.write(f"{i:6d} {a['name']:6s} {a['xyz'][0]:10.4f} {a['xyz'][1]:10.4f} "
                    f"{a['xyz'][2]:10.4f} {a['type']:4s} 1 LNK {a['charge']}\n")
        f.write("@<TRIPOS>BOND\n")
        for i, (a1, a2, bt) in enumerate(out_bonds, 1):
            f.write(f"{i:6d} {a1:6d} {a2:6d} {bt}\n")
        f.write("@<TRIPOS>SUBSTRUCTURE\n1 LNK 1 TEMP 0 **** **** 0 ROOT\n")
    print(f"{v}: {len(bcc_atoms)} -> {len(out_atoms)} atoms ({len(keep_heavy)} heavy), "
          f"{len(out_bonds)} bonds, netq={netq:+.3f}, NZ->C11={nz_c:.2f} A")
    from collections import Counter
    print(f"   types: {dict(Counter(a['type'] for a in out_atoms))}")

if __name__ == "__main__":
    for v in VARIANTS:
        rebuild(v)
