#!/usr/bin/env python3
"""Build clash-free idealized LNK initial geometries for exp-D (2026-07-17).

The MMFF conformers in the original bcc mol2 files are folded back through
the ECD surface (heavy-atom overlaps down to 0.5 A); rigid rotation about
the CE-NZ axis cannot untangle them. No RDKit available to regenerate
conformers, so the LNK coordinates are rebuilt with idealized chain
geometry:

- C11 anchored at NZ + 1.38 A along the CE->NZ axis (amide bond, unchanged)
- chain grown with ideal bond lengths/angles over the BFS tree; the main
  continuation (largest subtree) is placed first with a greedy dihedral
  scan (trans first, avoiding protein contacts <2.2 A, self-contacts
  <2.4 A, over-extension >22 A from NZ to keep the ~36k-atom box)
- sp2 centers stay planar (continuation trans; branch O trans to it)
- sp3 branch atoms locked tetrahedral at +/-120 deg from the continuation
- refinement sweeps re-scan continuations against the full current
  geometry until no 1-4+ self pair < 2.15 A remains
- H atoms added with standard tetrahedral/planar rules

Connectivity, atom names, GAFF2 types, charges and the NZ-C11 anchor are
identical to rebuild_lnk_mol2.py output -- only coordinates change.
Overwrites tleap/lnk_{v}_clean.mol2.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TLEAP = os.path.normpath(os.path.join(HERE, "..", "tleap"))
VARIANTS = ["no_linker", "gglu_1oeg", "gglu_2oeg", "gglu_3oeg"]

NZ_POS = np.array([-7.571, -0.448, 1.906])
CE_POS = np.array([-6.269, -0.196, 1.405])
AXIS = (NZ_POS - CE_POS) / np.linalg.norm(NZ_POS - CE_POS)

BOND_LEN = {frozenset(["c3", "c3"]): 1.54, frozenset(["c3", "c"]): 1.51,
            frozenset(["c", "o"]): 1.24, frozenset(["c", "n"]): 1.33,
            frozenset(["c3", "n"]): 1.46, frozenset(["c3", "n3"]): 1.47,
            frozenset(["c3", "os"]): 1.43}
ANGLE_AT = {"c3": 109.5, "n3": 109.5, "c": 120.0, "n": 120.0, "os": 111.0}
CANDIDATES = [180, 60, -60, 120, -120, 90, -90, 150, -150, 30, -30, 0]

def parse_mol2(path):
    atoms, bonds, sec = [], [], None
    for line in open(path):
        s = line.strip()
        if s.startswith("@<TRIPOS>"):
            sec = s; continue
        if sec == "@<TRIPOS>ATOM" and s:
            p = s.split()
            atoms.append(dict(name=p[1], type=p[5], charge=p[8]))
        elif sec == "@<TRIPOS>BOND" and s:
            p = s.split(); bonds.append((int(p[1]) - 1, int(p[2]) - 1, p[3]))
    return atoms, bonds

def parse_pdb_heavy(path):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                     for l in open(path)
                     if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() != "H"])

def nerf(A, B, C, r, theta_deg, phi_deg):
    """Place D given A-B-C, |C-D|=r, angle(B,C,D)=theta, dihedral(A,B,C,D)=phi."""
    th, ph = np.deg2rad(theta_deg), np.deg2rad(phi_deg)
    bc = C - B; bc /= np.linalg.norm(bc)
    ba = A - B
    n = np.cross(ba, bc)
    if np.linalg.norm(n) < 1e-9:
        tmp = np.array([1.0, 0, 0]) if abs(bc[0]) < 0.9 else np.array([0, 1.0, 0])
        n = np.cross(tmp, bc)
    n /= np.linalg.norm(n)
    n2 = np.cross(bc, n)
    return C + r * (-np.cos(th) * bc + np.sin(th) * (np.cos(ph) * n2 + np.sin(ph) * n))

def build(v, prot):
    path = os.path.join(TLEAP, f"lnk_{v}_clean.mol2")
    atoms, bonds = parse_mol2(path)
    n = len(atoms)
    heavy = [i for i, a in enumerate(atoms) if a["type"] not in ("hc", "hn")]
    hset = set(heavy)
    htype = {i: atoms[i]["type"] for i in heavy}
    adj = {i: [] for i in range(n)}
    for a1, a2, _ in bonds:
        adj[a1].append(a2); adj[a2].append(a1)

    def bond_len(t1, t2):
        return BOND_LEN.get(frozenset([t1, t2]), 1.45)

    c11 = next(i for i in heavy if atoms[i]["name"] == "C11")

    # BFS tree over heavy atoms, rooted at C11; children sorted by subtree size
    tree_par, tree_ch = {c11: None}, {i: [] for i in heavy}
    q = [c11]
    while q:
        cur = q.pop(0)
        for j in adj[cur]:
            if j in hset and j not in tree_par:
                tree_par[j] = cur; tree_ch[cur].append(j); q.append(j)
    def subtree(i):
        return 1 + sum(subtree(j) for j in tree_ch[i])
    for i in heavy:
        tree_ch[i].sort(key=lambda j: -subtree(j))

    pos = {c11: NZ_POS + 1.38 * AXIS}
    meta = {}      # j -> placement spec
    phi_cur = {}   # greedy continuation dihedrals
    order = []

    def apos(ref):
        return CE_POS if ref == "CE" else (NZ_POS if ref == "NZ" else pos[ref])

    def place(j, phi=None):
        m = meta[j]
        if m["kind"] == "greedy":
            phi = phi_cur[j] if phi is None else phi
        elif m["kind"] == "tetra":
            phi = (phi_cur[m["ref"]] + m["offset"]) % 360
        else:  # fixed absolute phi
            phi = m["phi"]
        pos[j] = nerf(apos(m["A"]), apos(m["B"]), pos[m["C"]], m["r"], m["ang"], phi)
        return phi

    def score_against(x, j):
        s = 0.0
        dp = np.sqrt(((prot - x) ** 2).sum(-1)).min()
        if dp < 2.2: s += 1000 * (2.2 - dp)
        s += -1.0 * min(dp, 6.0)  # soft clearance: steer the chain into solvent
        excl = {tree_par[j]} | set(adj[tree_par[j]]) | {j}
        for k, pk in pos.items():
            if k in excl or k not in hset: continue
            d = np.linalg.norm(pk - x)
            if d < 2.4: s += 500 * (2.4 - d)
        s += 0.5 * max(0.0, np.linalg.norm(x - NZ_POS) - 22.0)
        return s

    def frame(j, cur):
        """(A, B) refs for placing j about parent cur."""
        if cur == c11:
            return "CE", "NZ"
        pp = tree_par[cur]
        gp = tree_par[pp] if pp is not None else None
        return (gp if gp is not None else "CE"), pp

    # ---- initial placement over the BFS tree ----
    q = [c11]
    while q:
        cur = q.pop(0)
        cont_done = None
        for j in tree_ch[cur]:
            tc, tj = htype[cur], htype[j]
            r = bond_len(tc, tj)
            ang = ANGLE_AT.get(tc, 109.5)
            if cur == c11:
                ang = 123.0 if tj == "o" else 116.0  # frcmod N3-c-o / N3-c-c3
            A, B = frame(j, cur)
            if tj == "o" and tc == "c" and cont_done is None and htype[j] == "o" and not any(
                    htype[k] != "o" for k in tree_ch[cur]):
                # terminal COO- carbon: two O's, absolute phi 0 / 180
                n_o = sum(1 for k in tree_ch[cur] if k in meta)
                meta[j] = dict(kind="fixed", A=A, B=B, C=cur, r=r, ang=ang, phi=180.0 * n_o)
                place(j)
            elif tj == "o" and tc == "c":
                # sp2 branch O: trans to the continuation
                meta[j] = dict(kind="tetra", A=A, B=B, C=cur, r=r, ang=ang,
                               ref=cont_done, offset=180.0)
                place(j)
            elif cont_done is not None:
                # sp3 branch: tetrahedral +/-120 from the continuation
                offs = []
                for off in (120.0, -120.0):
                    meta[j] = dict(kind="tetra", A=A, B=B, C=cur, r=r, ang=ang,
                                   ref=cont_done, offset=off)
                    offs.append((score_against(place(j, 0), j), off))
                meta[j]["offset"] = min(offs)[1]
                place(j)
            else:
                # chain continuation (C13 is also greedy: the NZ-C11 amide
                # torsion relaxes to planar in the first ps via the
                # C8-N3-c-c3 dihedral; a forced trans start can trap the
                # chain inside a protein pocket)
                cand = CANDIDATES if cur == c11 else ([180] if tc in ("c", "n") else CANDIDATES)
                meta[j] = dict(kind="greedy", A=A, B=B, C=cur, r=r, ang=ang, cand=cand)
                best = min(cand, key=lambda p: score_against(
                    nerf(apos(A), apos(B), pos[cur], r, ang, p), j))
                phi_cur[j] = best
                place(j)
                cont_done = j
            order.append(j); q.append(j)

    # ---- refinement sweeps ----
    def worst_self():
        w = 99.0
        for ii, i in enumerate(heavy):
            for j in heavy[ii + 1:]:
                if j in adj[i] or any(j in adj[k] for k in adj[i]): continue
                w = min(w, np.linalg.norm(pos[i] - pos[j]))
        return w

    for sweep in range(20):
        if worst_self() >= 2.15: break
        for j in order:
            if meta[j]["kind"] != "greedy": continue
            m = meta[j]
            del pos[j]
            best = min(m["cand"], key=lambda p: score_against(
                nerf(apos(m["A"]), apos(m["B"]), pos[m["C"]], m["r"], m["ang"], p), j))
            phi_cur[j] = best
            place(j)
            for k2 in tree_ch[j]:
                place(k2)  # dependents follow (tetra/fixed)
        # safety: refresh everything downstream too
        for j in order:
            place(j)

    # ---- H placement ----
    for i, a in enumerate(atoms):
        if a["type"] not in ("hc", "hn"): continue
        hpars = [j for j in adj[i] if j in hset]
        assert len(hpars) == 1
        P = pos[hpars[0]]; pt = htype[hpars[0]]
        siblings = [k for k in adj[hpars[0]] if atoms[k]["type"] in ("hc", "hn")]
        ksis = siblings.index(i)
        heavy_nbrs = [pos[k] for k in adj[hpars[0]] if k in hset and k != i]
        if pt == "c3" and len(heavy_nbrs) == 2:
            u1 = heavy_nbrs[0] - P; u1 /= np.linalg.norm(u1)
            u2 = heavy_nbrs[1] - P; u2 /= np.linalg.norm(u2)
            p = np.cross(u1, u2); p /= np.linalg.norm(p)
            d = -0.5 * (u1 + u2) + (np.sqrt(2 / 3) if ksis == 0 else -np.sqrt(2 / 3)) * p
            pos[i] = P + 1.09 * d / np.linalg.norm(d)
        elif pt == "c3" and len(heavy_nbrs) == 3:
            d = -sum((x - P) / np.linalg.norm(x - P) for x in heavy_nbrs)
            pos[i] = P + 1.09 * d / np.linalg.norm(d)
        elif pt in ("n", "n3") and len(heavy_nbrs) == 2:
            d = -sum((x - P) / np.linalg.norm(x - P) for x in heavy_nbrs)
            pos[i] = P + 1.01 * d / np.linalg.norm(d)
        elif pt in ("n", "n3") and len(heavy_nbrs) == 1:
            u = (heavy_nbrs[0] - P) / np.linalg.norm(heavy_nbrs[0] - P)
            tmp = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
            p = np.cross(u, tmp); p /= np.linalg.norm(p)
            q = np.cross(u, p)
            az = np.deg2rad(90 if ksis == 0 else 270)
            d = np.cos(np.deg2rad(109.5)) * u + np.sin(np.deg2rad(109.5)) * (np.cos(az) * p + np.sin(az) * q)
            pos[i] = P + 1.01 * d / np.linalg.norm(d)
        else:
            raise RuntimeError(f"{v}: unhandled H parent {atoms[hpars[0]]['name']} ({pt})")

    # ---- validation ----
    xyz = np.array([pos[i] for i in range(n)])
    bad = [(np.linalg.norm(xyz[a1] - xyz[a2]), atoms[a1]["name"], atoms[a2]["name"])
           for a1, a2, _ in bonds if not (0.9 < np.linalg.norm(xyz[a1] - xyz[a2]) < 1.7)]
    assert not bad, f"{v}: bad bond lengths {bad[:5]}"
    w_self, w_pair = 99.0, None
    for ii, i in enumerate(heavy):
        for j in heavy[ii + 1:]:
            if j in adj[i] or any(j in adj[k] for k in adj[i]): continue
            d = np.linalg.norm(xyz[i] - xyz[j])
            if d < w_self: w_self, w_pair = d, (atoms[i]["name"], atoms[j]["name"])
    dmat = np.sqrt(((prot[None, :, :] - xyz[heavy][:, None, :]) ** 2).sum(-1))
    dmin, nclash = 99.0, 0
    for ii in range(len(heavy)):
        for k in range(len(prot)):
            if np.linalg.norm(prot[k] - NZ_POS) < 0.01: continue  # the NZ-C11 bond
            d = dmat[ii, k]
            if d < 1.9: nclash += 1
            dmin = min(dmin, d)
    print(f"{v}: NZ-C11={np.linalg.norm(pos[c11]-NZ_POS):.2f} A, prot min={dmin:.2f} A "
          f"(contacts<1.9A={nclash}), self min(1-4+)={w_self:.2f} A ({w_pair[0]}-{w_pair[1]}), "
          f"max r(NZ)={max(np.linalg.norm(pos[i]-NZ_POS) for i in heavy):.1f} A")
    assert nclash == 0 and dmin > 1.9 and w_self >= 2.1, f"{v}: geometry validation failed"

    with open(path, "w") as f:
        f.write("@<TRIPOS>MOLECULE\nLNK\n")
        f.write(f" {n:5d} {len(bonds):5d}     0     0     0\nSMALL\nGAFF2\n")
        f.write(f"rebuilt {v} linker-C18 (idealized clash-free geometry)\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for i, a in enumerate(atoms, 1):
            f.write(f"{i:6d} {a['name']:6s} {xyz[i-1][0]:10.4f} {xyz[i-1][1]:10.4f} "
                    f"{xyz[i-1][2]:10.4f} {a['type']:4s} 1 LNK {a['charge']}\n")
        f.write("@<TRIPOS>BOND\n")
        for i, (a1, a2, bt) in enumerate(bonds, 1):
            f.write(f"{i:6d} {a1 + 1:6d} {a2 + 1:6d} {bt}\n")
        f.write("@<TRIPOS>SUBSTRUCTURE\n1 LNK 1 TEMP 0 **** **** 0 ROOT\n")

if __name__ == "__main__":
    prot = parse_pdb_heavy(os.path.join(TLEAP, "ecd_pep_nocap.pdb"))
    for v in VARIANTS:
        build(v, prot)
    print("Done.")
