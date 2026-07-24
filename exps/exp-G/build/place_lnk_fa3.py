#!/usr/bin/env python3
"""Place LNK (gGlu-2xOEG-C18 diacid) with its tail pre-anchored at HSA FA3.

Anchoring geometry (exp-G design, see README.md):
- Distal carboxyl C55/O56/O57 on the 1E7G MYR 1003 carboxyl (crystal plane)
  -> double salt bridge ARG346+ARG483 (crystal numbering ARG348/485)
- C54..C45 along crystal MYR C2..C11 through the pocket
- C44..C37 along the exp-C validated FA3 exit path (c18_diacid_fa3.mol2
  ext1..ext7, all clearances >= 2.2 A); O38 dihedral scanned toward ARG410
  (crystal numbering = topology ARG408)
- ADO2/ADO1/gGlu grown greedily into solvent (hard clash penalty + local
  exit direction away from the protein surface)
- O12 dihedral chosen so the implied NZ direction points into open solvent
- NZ target = C11 + 1.38 A, planar amide (trans to O12)

Coordinates are built fresh; atom names/types/charges/bonds are inherited
verbatim from exps/exp-D/tleap/lnk_gglu_2oeg_clean.mol2.
Output: exps/exp-G/tleap/lnk_2oeg_fa3.mol2 + nz_target.txt

Frame convention: nerf(A, B, C, r, ang, phi) places a new atom bonded to C
with angle(B,C,new)=ang and dihedral(A,B,C,new)=phi. For chain atom
chain[k], the frame is (chain[k-3], chain[k-2], chain[k-1]); for a branch
O on carbonyl chain[k], the frame is that of the continuation chain[k+1],
i.e. (chain[k-2], chain[k-1], chain[k]).
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXP_G = os.path.normpath(os.path.join(HERE, ".."))
TLEAP_D = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap"
PDB_1E7G = "/home/scroll/personal/semaglutide-in-silico/exps/exp-C/structures/1E7G.pdb"
HSA_PDB = "/home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap/hsa_no_myr.pdb"
EXPC_MOL2 = "/home/scroll/personal/semaglutide-in-silico/exps/exp-C/tleap/c18_diacid_fa3.mol2"
OUT = os.path.join(EXP_G, "tleap", "lnk_2oeg_fa3.mol2")

BOND_LEN = {frozenset(["c3", "c3"]): 1.54, frozenset(["c3", "c"]): 1.51,
            frozenset(["c", "o"]): 1.24, frozenset(["c", "n"]): 1.33,
            frozenset(["c", "c3"]): 1.51,
            frozenset(["c3", "n"]): 1.46, frozenset(["c3", "n3"]): 1.47,
            frozenset(["c3", "os"]): 1.43}
ANGLE_AT = {"c3": 109.5, "n3": 109.5, "c": 120.0, "n": 120.0, "os": 111.0}
CANDIDATES = [180, 60, -60, 120, -120, 90, -90, 150, -150, 30, -30, 0]

def parse_mol2(path, with_coords=False):
    atoms, bonds, sec = [], [], None
    for line in open(path):
        s = line.strip()
        if s.startswith("@<TRIPOS>"):
            sec = s; continue
        if sec == "@<TRIPOS>ATOM" and s:
            p = s.split()
            a = dict(id=int(p[0]), name=p[1], type=p[5], charge=p[8])
            if with_coords:
                a["xyz"] = np.array([float(p[2]), float(p[3]), float(p[4])])
            atoms.append(a)
        elif sec == "@<TRIPOS>BOND" and s:
            p = s.split(); bonds.append((int(p[1]) - 1, int(p[2]) - 1, p[3]))
    return (atoms, bonds) if not with_coords else atoms

def load_myr():
    atoms = {}
    for line in open(PDB_1E7G):
        if line.startswith("HETATM") and line[17:20].strip() == "MYR" and int(line[22:26]) == 1003:
            atoms[line[12:16].strip()] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return atoms

def load_hsa_heavy():
    pts, res = [], {}
    for line in open(HSA_PDB):
        if line.startswith("ATOM"):
            name = line[12:16].strip(); rn = int(line[22:26])
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            pts.append(xyz); res.setdefault(rn, {})[name] = xyz
    return np.array(pts), res

def nerf(A, B, C, r, theta_deg, phi_deg):
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

def unit(v):
    return v / np.linalg.norm(v)

MYR = load_myr()
HSA, HSA_RES = load_hsa_heavy()
ARG410 = HSA_RES[410]  # crystal numbering = topology ARG408
ARG410_NH = np.array([ARG410["NH1"], ARG410["NH2"]])
EXPC = {a["name"]: a["xyz"] for a in parse_mol2(EXPC_MOL2, with_coords=True)}

atoms, bonds = parse_mol2(os.path.join(TLEAP_D, "lnk_gglu_2oeg_clean.mol2"))
N = len(atoms)
name2i = {a["name"]: i for i, a in enumerate(atoms)}
htype = {i: atoms[i]["type"] for i in range(N)}
adj = {i: set() for i in range(N)}
for a1, a2, _ in bonds:
    adj[a1].add(a2); adj[a2].add(a1)
heavy = [i for i in range(N) if htype[i] not in ("hc", "hn")]
hset = set(heavy)

def bfs_path(s, t):
    prev = {s: None}; q = [s]
    while q:
        u = q.pop(0)
        if u == t: break
        for w in sorted(adj[u]):
            if w in hset and w not in prev:
                prev[w] = u; q.append(w)
    path = [t]
    while path[-1] != s:
        path.append(prev[path[-1]])
    return path[::-1]

chain = bfs_path(name2i["C55"], name2i["C11"])
assert len(chain) == 40
print(f"main chain: {'-'.join(atoms[i]['name'] for i in chain)}")

pos = {}

# ---- 1. distal carboxyl on the MYR carboxyl (crystal plane) ----
pos[name2i["C55"]] = MYR["C1"].copy()
u = unit(MYR["C2"] - MYR["C1"])
w = MYR["O1"] - MYR["C1"]; w -= np.dot(w, u) * u; w = unit(w)
for oname, sign in (("O56", +1), ("O57", -1)):
    d = np.cos(np.deg2rad(118)) * u + sign * np.sin(np.deg2rad(118)) * w
    pos[name2i[oname]] = MYR["C1"] + 1.26 * d
print(f"O56->MYR O1: {np.linalg.norm(pos[name2i['O56']]-MYR['O1']):.2f} A, "
      f"O57->MYR O2: {np.linalg.norm(pos[name2i['O57']]-MYR['O2']):.2f} A")

# ---- 2. C54..C45 = MYR C2..C11 ----
for k in range(10):
    pos[name2i[f"C{54-k}"]] = MYR[f"C{2+k}"].copy()

# ---- 3. C44..C39 + C37 = exp-C validated exit path ext1..ext7 (C11..C17) ----
for nm, src in (("C44", "C11"), ("C43", "C12"), ("C42", "C13"),
                ("C41", "C14"), ("C40", "C15"), ("C39", "C16"), ("C37", "C17")):
    pos[name2i[nm]] = EXPC[src].copy()
    d = np.sqrt(((HSA - pos[name2i[nm]]) ** 2).sum(-1)).min()
    assert d >= 2.2, f"{nm}<-expC {src}: clearance {d:.2f}"
# bond sanity along the mapped segment
for a, b in (("C45", "C44"), ("C44", "C43"), ("C43", "C42"), ("C42", "C41"),
             ("C41", "C40"), ("C40", "C39"), ("C39", "C37")):
    d = np.linalg.norm(pos[name2i[a]] - pos[name2i[b]])
    assert 1.3 < d < 1.8, f"{a}-{b} = {d:.2f}"
print("C44..C37 mapped from exp-C exit path (clearances >= 2.2 A)")

# ---- greedy machinery ----
EXIT_LOCAL = [None]  # set after C36; guide direction for the ADO segment

def score(x, excl):
    s = 0.0
    dp = np.sqrt(((HSA - x) ** 2).sum(-1)).min()
    if dp < 2.6: s += 2000 * (2.6 - dp)
    s += -1.5 * min(dp, 8.0)
    if EXIT_LOCAL[0] is not None:
        s += -0.3 * min(np.dot(x - pos[name2i["C37"]], EXIT_LOCAL[0]), 25.0)
    for j in heavy:
        if j in excl or j not in pos: continue
        d = np.linalg.norm(pos[j] - x)
        if d < 2.4: s += 500 * (2.4 - d)
    return s

def grow(k, cand=CANDIDATES, extra=None):
    """Place chain[k] using frame (chain[k-3], chain[k-2], chain[k-1])."""
    i = chain[k]
    par, gp, ggp = chain[k - 1], chain[k - 2], chain[k - 3]
    r = BOND_LEN.get(frozenset([htype[par], htype[i]]), 1.45)
    ang = ANGLE_AT.get(htype[par], 109.5)
    best, best_phi = None, None
    for phi in cand:
        x = nerf(pos[ggp], pos[gp], pos[par], r, ang, phi)
        sc = score(x, {i, par, gp, ggp} | set(adj[par]))
        if extra: sc += extra(x, phi)
        if best is None or sc < best:
            best, best_phi = sc, phi
    pos[i] = nerf(pos[ggp], pos[gp], pos[par], r, ang, best_phi)
    return best_phi

def place_branch_o(k_c, i_o, phi_cont):
    """sp2 branch O on carbonyl chain[k_c], trans to continuation chain[k_c+1]."""
    iC = chain[k_c]
    r = BOND_LEN[frozenset(["c", "o"])]
    pos[i_o] = nerf(pos[chain[k_c - 2]], pos[chain[k_c - 1]], pos[iC], r, 120.0,
                    (phi_cont + 180) % 360)

# ---- 4. O38 toward ARG410; C36 planar ----
k37 = 17
i37 = chain[k37]
best = None
for phi in range(0, 360, 5):
    x = nerf(pos[chain[k37 - 2]], pos[chain[k37 - 1]], pos[i37], 1.24, 120.0, phi)
    dnh = np.sqrt(((ARG410_NH - x) ** 2).sum(-1)).min()
    sc = abs(dnh - 3.0) * 2.0 + (2000 * (2.6 - np.sqrt(((HSA - x) ** 2).sum(-1)).min())
                                 if np.sqrt(((HSA - x) ** 2).sum(-1)).min() < 2.6 else 0)
    if best is None or sc < best:
        best, phi_o38 = sc, phi
pos[name2i["O38"]] = nerf(pos[chain[k37 - 2]], pos[chain[k37 - 1]], pos[i37], 1.24, 120.0, phi_o38)
d_o38 = np.sqrt(((ARG410_NH - pos[name2i["O38"]]) ** 2).sum(-1)).min()
print(f"C37->ARG410 CZ: {np.linalg.norm(pos[i37]-ARG410['CZ']):.2f} A, O38->ARG410 NHx: {d_o38:.2f} A")
# C36 (idx 18): same frame as O38, trans
pos[chain[18]] = nerf(pos[chain[k37 - 2]], pos[chain[k37 - 1]], pos[i37],
                      BOND_LEN[frozenset(["c", "c3"])], 120.0, (phi_o38 + 180) % 360)

# local exit direction: away from the protein centroid near C37
near = HSA[np.sqrt(((HSA - pos[i37]) ** 2).sum(-1)) < 15.0]
EXIT_LOCAL[0] = unit(pos[i37] - near.mean(axis=0))
print(f"exit direction: {np.round(EXIT_LOCAL[0], 2)} (from {len(near)}-atom local centroid)")

# ---- 5. chain idx 19..38 (O35..C13) ----
phi_used = {}
for k in range(19, 39):
    nm = atoms[chain[k]]["name"]
    if nm == "C26":
        phi_used[k] = grow(k, cand=[180])
        place_branch_o(k - 1, name2i["O28"], phi_used[k])
    elif nm == "C15":
        phi_used[k] = grow(k, cand=[180])
        place_branch_o(k - 1, name2i["O18"], phi_used[k])
    else:
        phi_used[k] = grow(k)

# N16 (gGlu NH2): tetrahedral branch on C15, offset +/-120 from C14's dihedral
k14 = next(k for k in range(len(chain)) if atoms[chain[k]]["name"] == "C14")
i15, i16 = name2i["C15"], name2i["N16"]
offs = []
for off in (120.0, -120.0):
    x = nerf(pos[chain[k14 - 3]], pos[chain[k14 - 2]], pos[i15],
             BOND_LEN[frozenset(["c3", "n3"])], 109.5, (phi_used[k14] + off) % 360)
    offs.append((score(x, {i16, i15, chain[k14-2], chain[k14-3]} | set(adj[i15])), off))
best_off = min(offs)[1]
pos[i16] = nerf(pos[chain[k14 - 3]], pos[chain[k14 - 2]], pos[i15],
                BOND_LEN[frozenset(["c3", "n3"])], 109.5, (phi_used[k14] + best_off) % 360)

# ---- 6. C11 (idx 39) greedy; joint (NZ, ECD-rotation) docking search ----
k11 = 39
grow(k11)
i11 = chain[k11]

# ECD-peptide template (rigid 3IOL pose): heavy coords + Lys26 anchors
ECD_PEP = "/home/scroll/personal/semaglutide-in-silico/exps/exp-D/tleap/ecd_pep_nocap.pdb"
ecd_meta, ecd_xyz = [], []
for line in open(ECD_PEP):
    if line.startswith("ATOM"):
        ecd_meta.append((line[21], int(line[22:26]), line[12:16].strip()))
        ecd_xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
ecd_xyz = np.array(ecd_xyz)
nz0 = ce0 = None
for i, (ch, rn, nm) in enumerate(ecd_meta):
    if ch == "B" and rn == 117 and nm == "NZ": nz0 = ecd_xyz[i]
    if ch == "B" and rn == 117 and nm == "CE": ce0 = ecd_xyz[i]
assert nz0 is not None and ce0 is not None
lys26_idx = {i for i, (ch, rn, nm) in enumerate(ecd_meta) if ch == "B" and rn == 117}
ecd_core_idx = [i for i, (ch, rn, nm) in enumerate(ecd_meta) if ch == "A"]

def rot_align(v0, v1):
    k = np.cross(v0, v1)
    if np.linalg.norm(k) < 1e-12:
        return np.eye(3) if np.dot(v0, v1) > 0 else -np.eye(3)
    k /= np.linalg.norm(k)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    ang = np.arccos(np.clip(np.dot(v0, v1), -1, 1))
    return np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K

def rot_axis(axis, theta):
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * K @ K

c11p = pos[name2i["C11"]]
lnk_placed = np.array([pos[j] for j in heavy if j in pos], dtype=float)
obstacles = np.vstack([HSA, lnk_placed])
hsa_cent = HSA.mean(axis=0)

best = None
for phi in range(0, 360, 5):
    nz = nerf(pos[chain[k11 - 2]], pos[chain[k11 - 1]], pos[i11], 1.38, 116.0, (phi + 180) % 360)
    axis = unit(nz - c11p)
    R0 = rot_align(unit(ce0 - nz0), axis)  # NZ->CE aligns with C11->NZ (CE beyond NZ)
    e1 = (R0 @ (ecd_xyz - nz0).T).T + nz
    cent0 = e1[ecd_core_idx].mean(axis=0)
    for deg in range(0, 360, 5):
        R = rot_axis(axis, np.deg2rad(deg))
        xyz = (R @ (e1 - nz).T).T + nz
        d = np.sqrt(((obstacles[None, :, :] - xyz[:, None, :]) ** 2).sum(-1))
        nc = 0
        for ii in range(len(xyz)):
            if ii in lys26_idx: continue
            nc += int((d[ii] < 2.2).sum())
        dmin = d.min()
        sep = np.linalg.norm(((R @ (cent0 - nz).T).T + nz) - hsa_cent)
        key = (nc == 0, round(dmin, 2), round(sep, 1))
        if best is None or key > best[0]:
            best = (key, phi, deg, nz.copy(), axis.copy(), R0.copy())
(ok, dmin, sep), phi_o12, deg, nz_target, axis, R0 = best
print(f"joint docking search: phi_O12={phi_o12}, rot={deg}, clash-free={ok}, "
      f"min dist {dmin:.2f} A, ECD-HSA sep {sep:.1f} A")
assert ok, "no (NZ, rotation) pair eliminates clashes"

pos[name2i["O12"]] = nerf(pos[chain[k11 - 2]], pos[chain[k11 - 1]], pos[i11], 1.24, 120.0, phi_o12)
print(f"NZ target: {np.round(nz_target, 2)}, span C55->NZ: "
      f"{np.linalg.norm(nz_target - pos[name2i['C55']]):.1f} A")
np.savez(os.path.join(EXP_G, "tleap", "ecd_dock_transform.npz"),
         nz_target=nz_target, axis=axis, R0=R0, deg=float(deg))

# ---- 7. H placement (exp-D rules) ----
for i, a in enumerate(atoms):
    if htype[i] not in ("hc", "hn"): continue
    par = next(j for j in adj[i] if j in hset)
    P = pos[par]; pt = htype[par]
    sib = [k for k in adj[par] if htype[k] in ("hc", "hn")]
    ks = sib.index(i)
    nbrs = [pos[k] for k in adj[par] if k in hset]
    if pt == "c3" and len(nbrs) == 2:
        u1 = unit(nbrs[0] - P); u2 = unit(nbrs[1] - P)
        p = unit(np.cross(u1, u2))
        d = -0.5 * (u1 + u2) + (np.sqrt(2 / 3) if ks == 0 else -np.sqrt(2 / 3)) * p
        pos[i] = P + 1.09 * unit(d)
    elif pt == "c3" and len(nbrs) == 3:
        d = -sum(unit(x - P) for x in nbrs)
        pos[i] = P + 1.09 * unit(d)
    elif pt in ("n", "n3") and len(nbrs) == 2:
        d = -sum(unit(x - P) for x in nbrs)
        pos[i] = P + 1.01 * unit(d)
    elif pt in ("n", "n3") and len(nbrs) == 1:
        u = unit(nbrs[0] - P)
        tmp = np.array([1.0, 0, 0]) if abs(u[0]) < 0.9 else np.array([0, 1.0, 0])
        p = unit(np.cross(u, tmp)); q = np.cross(u, p)
        az = np.deg2rad(90 if ks == 0 else 270)
        d = np.cos(np.deg2rad(109.5)) * u + np.sin(np.deg2rad(109.5)) * (np.cos(az) * p + np.sin(az) * q)
        pos[i] = P + 1.01 * unit(d)

# ---- 8. validation + write ----
xyz = np.array([pos[i] for i in range(N)])
bad = [(np.linalg.norm(xyz[a1] - xyz[a2]), atoms[a1]["name"], atoms[a2]["name"])
       for a1, a2, _ in bonds if not (0.9 < np.linalg.norm(xyz[a1] - xyz[a2]) < 1.7)]
assert not bad, f"bad bond lengths {bad[:5]}"
w_self, w_pair = 99.0, None
for ii, i in enumerate(heavy):
    for j in heavy[ii + 1:]:
        if j in adj[i] or any(j in adj[k] for k in adj[i]): continue
        d = np.linalg.norm(xyz[i] - xyz[j])
        if d < w_self: w_self, w_pair = d, (atoms[i]["name"], atoms[j]["name"])
d_hsa = np.sqrt(((HSA[None, :, :] - xyz[heavy][:, None, :]) ** 2).sum(-1)).min()
print(f"self min(1-4+): {w_self:.2f} A ({w_pair[0]}-{w_pair[1]}), LNK-HSA min: {d_hsa:.2f} A")
assert w_self >= 2.0 and d_hsa > 1.9

os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(OUT, "w") as f:
    f.write("@<TRIPOS>MOLECULE\nLNK\n")
    f.write(f" {N:5d} {len(bonds):5d}     0     0     0\nSMALL\nGAFF2\n")
    f.write("gGlu-2xOEG-C18 pre-anchored at HSA FA3 (exp-G)\n\n@<TRIPOS>ATOM\n")
    for i, a in enumerate(atoms, 1):
        f.write(f"{i:6d} {a['name']:6s} {xyz[i-1][0]:10.4f} {xyz[i-1][1]:10.4f} "
                f"{xyz[i-1][2]:10.4f} {a['type']:4s} 1 LNK {a['charge']}\n")
    f.write("@<TRIPOS>BOND\n")
    for i, (a1, a2, bt) in enumerate(bonds, 1):
        f.write(f"{i:6d} {a1 + 1:6d} {a2 + 1:6d} {bt}\n")
    f.write("@<TRIPOS>SUBSTRUCTURE\n1 LNK 1 TEMP 0 **** **** 0 ROOT\n")
np.savetxt(os.path.join(EXP_G, "tleap", "nz_target.txt"), nz_target, fmt="%.4f")
print(f"wrote {OUT} and nz_target.txt")
