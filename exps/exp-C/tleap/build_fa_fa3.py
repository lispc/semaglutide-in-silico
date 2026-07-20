#!/usr/bin/env python3
"""Build all-atom GAFF2 mol2 for free fatty acids positioned at the HSA FA3 site.

Supports arbitrary chain length (v4, 2026-07-20):
  Usage: python build_fa_fa3.py {mono|diacid|both} L [L ...] [--force]
    L = TOTAL carbon count (chemical convention):
      diacid L -> HOOC-(CH2){L-2}-COOH  (e.g. C18 diacid = 16 CH2, 54 atoms)
      mono   L -> CH3-(CH2){L-2}-COO-   (e.g. C18 stearate = 17 chain C, 56 atoms)
    Output: c{L}_{monoacid|diacid}_fa3.mol2 ; existing files need --force.

  NOTE (convention flag): the 2026-05-27 "c18" controls were built with 18
  CHAIN carbons (C01..C18), i.e. mono = nonadecanoate (19 total C, 58 atoms),
  diacid = eicosanedioate (20 total C, 60 atoms) -- NOT true C18. Reproduce
  those molecules with:  mono 19  /  diacid 20 .

Placement strategy (v3, 2026-07-17):
  - C1X/C1D (distal carboxyl C) sits exactly on MYR 1003 C1 (1E7G FA3);
    O1D/O2D are built as ideal sp2 carboxyl in the MYR carboxyl plane,
    so they land on the crystal oxygens (~2.8 A salt bridge to ARG348/485).
  - Chain C01..C10 sit on the crystallographic MYR 1003 C2..C11 positions.
  - C11..Cn (and proximal C1P) extend the chain with a greedy dihedral
    scan per bond (full 360 deg, prefer near-trans) that maximizes the
    minimum distance to HSA heavy atoms -> clash-free curved tail, matching
    the crystal observation that MYR C12+ is disordered (roomy pocket end).
  - Shorter chains simply END earlier (tail buried shallower) -- the chain
    is never artificially stretched or folded.
  (The original builders rotated +112 deg about a FIXED axis per bond,
  giving a 68 deg internal angle and a self-coiled chain -- C01/C17
  overlapped at 0.13 A; and a naive straight zigzag clashes with PHE403.)

Charge scheme (consistent with the 2026-05-27 controls):
  GAFF2 types c/o/c3/hc, raw c=+0.70 o=-0.80 c3=-0.12 hc=+0.06 uniformly
  scaled to integer total (-1 mono / -2 diacid), rounded to 4 decimals with
  the residual absorbed on a mid-chain carbon.
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PDB = os.path.join(HERE, "..", "structures", "1E7G.pdb")
HSA_PDB = os.path.join(HERE, "hsa_no_myr.pdb")

BOND_CC, BOND_CH, BOND_CO = 1.54, 1.09, 1.26
ANG_INT = 112 * np.pi / 180          # internal C-C-C angle
ANG_TURN = np.pi - ANG_INT           # 68 deg between successive bond vectors
ANG_CCO = 118 * np.pi / 180          # carboxyl C-C-O angle

# ---------------- template loading ----------------
def load_myr():
    atoms = {}
    with open(PDB) as f:
        for line in f:
            if line.startswith("HETATM") and line[17:20].strip() == "MYR" \
               and int(line[22:26]) == 1003:
                atoms[line[12:16].strip()] = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])])
    need = ["C1", "O1", "O2"] + [f"C{i}" for i in range(2, 12)]
    missing = [n for n in need if n not in atoms]
    if missing:
        raise RuntimeError(f"MYR 1003 incomplete in {PDB}: missing {missing}")
    return atoms

def load_hsa_heavy():
    pts = []
    with open(HSA_PDB) as f:
        for line in f:
            if line.startswith("ATOM") and line[76:78].strip() != "H":
                pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return np.array(pts)

MYR = load_myr()
HSA = load_hsa_heavy()

def unit(v):
    return v / np.linalg.norm(v)

def rot_about(v, axis, ang):
    """Rodrigues rotation of v about unit axis."""
    return (v * np.cos(ang) + np.cross(axis, v) * np.sin(ang)
            + axis * np.dot(axis, v) * (1 - np.cos(ang)))

class Atom:
    def __init__(self, name, pos, atype, charge):
        self.name = name; self.pos = np.array(pos, dtype=float)
        self.atype = atype; self.charge = charge

# ---------------- geometry builders ----------------
def carboxyl(cname, onames, ccx, u_to_chain, w_hint):
    """Ideal sp2 carboxyl: C-O 1.26 A, C-C-O 118 deg (O-C-O 124 deg).
    w_hint selects the O-C-O plane direction (component perp to u_to_chain)."""
    w = w_hint - np.dot(w_hint, u_to_chain) * u_to_chain
    w = w / np.linalg.norm(w)
    d1 = np.cos(ANG_CCO) * u_to_chain + np.sin(ANG_CCO) * w
    d2 = np.cos(ANG_CCO) * u_to_chain - np.sin(ANG_CCO) * w
    return [Atom(cname, ccx, 'c', 0.70),
            Atom(onames[0], ccx + d1 * BOND_CO, 'o', -0.80),
            Atom(onames[1], ccx + d2 * BOND_CO, 'o', -0.80)]

def next_bond_dir(p3, p2, p1, phi):
    """Direction for the next bond p1->p_new such that the internal angle at
    p1 is 112 deg (turn 68 deg from unit(p2->p1)); phi rotates the new bond
    about the (p2->p1) axis, phi=0 = trans (extended planar zigzag)."""
    u_prev2 = unit(p2 - p3)      # bond p3->p2
    u_prev = unit(p1 - p2)       # bond p2->p1
    n = np.cross(u_prev2, u_prev)
    ln = np.linalg.norm(n)
    if ln < 1e-6:
        n = np.cross(u_prev, [0., 0., 1.])
        if np.linalg.norm(n) < 1e-6:
            n = np.cross(u_prev, [1., 0., 0.])
    n = n / np.linalg.norm(n)
    d = np.cos(ANG_TURN) * u_prev + np.sin(ANG_TURN) * np.cross(n, u_prev)
    # choose trans reference: new point farther from p3 (extended, not folded)
    d_alt = np.cos(ANG_TURN) * u_prev - np.sin(ANG_TURN) * np.cross(n, u_prev)
    if np.linalg.norm(p1 + d_alt * BOND_CC - p3) > np.linalg.norm(p1 + d * BOND_CC - p3):
        d = d_alt
    return rot_about(d, u_prev, phi)

def extend_chain(c_pos, n_ext, extra_pts, label):
    """Greedily append n_ext carbons; full-circle dihedral scan per bond,
    maximizing min distance to (HSA heavy + already-placed FA heavy)."""
    placed = [a.pos for a in extra_pts]
    for k in range(n_ext):
        p3, p2, p1 = c_pos[-3], c_pos[-2], c_pos[-1]
        best = None
        obstacles = c_pos[:-2] + placed  # exclude 1-2/1-3 neighbors
        obstacles = np.array(obstacles) if obstacles else None
        for deg in range(0, 360, 5):
            phi = np.deg2rad(deg)
            d = next_bond_dir(p3, p2, p1, phi)
            cand = p1 + d * BOND_CC
            dmin = np.linalg.norm(HSA - cand, axis=1).min()
            if obstacles is not None and len(obstacles):
                dmin = min(dmin, np.linalg.norm(obstacles - cand, axis=1).min())
            # prefer near-trans (phi ~ 0) on ties
            score = dmin - 0.05 * min(deg, 360 - deg) / 180.0
            if best is None or score > best[0]:
                best = (score, dmin, deg, cand)
        _, dmin, deg, cand = best
        if dmin < 1.5:
            print(f"  WARNING [{label}]: extension C at {deg} deg still only {dmin:.2f} A from environment")
        c_pos.append(cand)
        print(f"  [{label}] ext C{len(c_pos):02d}: dihedral offset {deg:3d} deg, min dist {dmin:.2f} A")
    return c_pos

def add_h(hydro, name, cp, v1, v2):
    """Two tetrahedral H on a CH2 carbon given its two C-neighbor vectors."""
    v1u, v2u = unit(v1), unit(v2)
    bis = -(v1u + v2u)
    bis = bis / np.linalg.norm(bis)
    perp_h = np.cross(v1u, v2u)
    if np.linalg.norm(perp_h) < 0.01:
        perp_h = np.array([0., 1., 0.])
    perp_h = perp_h / np.linalg.norm(perp_h)
    hydro.append(Atom(f'{name}A', cp + unit(bis + perp_h) * BOND_CH, 'hc', 0.06))
    hydro.append(Atom(f'{name}B', cp + unit(bis - perp_h) * BOND_CH, 'hc', 0.06))

# ---------------- main build ----------------
def build(kind, n_chain, label):
    """kind: 'mono' (charge -1) or 'diacid' (charge -2).
    n_chain: number of methylene-chain carbons C01..Cn.
    label: system label used for the output filename and reports."""
    if n_chain < 3:
        raise ValueError(f"n_chain={n_chain} too small")
    heavy, hydro = [], []

    # Distal carboxyl exactly on the MYR carboxyl (crystal plane via O1 hint)
    ccx_d = MYR["C1"]
    c01 = MYR["C2"]
    u_to_chain = unit(c01 - ccx_d)
    heavy += carboxyl('C1X' if kind == 'mono' else 'C1D', ('O1D', 'O2D'),
                      ccx_d, u_to_chain, MYR["O1"] - ccx_d)

    # C01..C10 on crystallographic MYR C2..C11 (or fewer for short chains)
    n_seed = min(n_chain, 10)
    c_pos = [MYR[f"C{i}"] for i in range(2, 2 + n_seed)]

    # Remaining chain carbons: greedy clash-free extension
    c_pos = extend_chain(list(c_pos), n_chain - n_seed, heavy[:3], label)

    for i, cp in enumerate(c_pos):
        heavy.append(Atom(f'C{i+1:02d}', cp, 'c3', -0.12))

    # Proximal carboxyl (diacid): clash-aware continuation from the last
    # chain carbon, scored over the WHOLE carboxyl group (C1P + O1P + O2P)
    if kind == 'diacid':
        heavy_arr = np.array([a.pos for a in heavy])
        best = None
        for deg in range(0, 360, 5):
            b_next = next_bond_dir(c_pos[-3], c_pos[-2], c_pos[-1], np.deg2rad(deg))
            ccx_p = c_pos[-1] + b_next * 1.52
            w = np.cross(b_next, unit(c_pos[-1] - c_pos[-2]))
            grp = carboxyl('C1P', ('O1P', 'O2P'), ccx_p, -b_next, w)
            pts = np.array([a.pos for a in grp])
            dmin = min(np.linalg.norm(HSA[None, :, :] - pts[:, None, :], axis=-1).min(),
                       np.linalg.norm(heavy_arr[None, :, :] - pts[:, None, :], axis=-1).min())
            score = dmin - 0.05 * min(deg, 360 - deg) / 180.0
            if best is None or score > best[0]:
                best = (score, dmin, deg, grp)
        _, dmin, deg, grp = best
        heavy += grp
        print(f"  [{label}] C1P carboxyl: offset {deg} deg, min dist {dmin:.2f} A")

    # Hydrogens
    for i in range(n_chain):
        cp = c_pos[i]
        if kind == 'mono' and i == n_chain - 1:
            # Terminal CH3: proper tetrahedral geometry (H-C-C = 109.5 deg)
            back = unit(c_pos[i-1] - cp)
            p2 = np.cross(back, [0., 0., 1.])
            if np.linalg.norm(p2) < 0.01:
                p2 = np.cross(back, [1., 0., 0.])
            p2 = p2 / np.linalg.norm(p2)
            p1 = np.cross(p2, back)
            base = cp - back * (BOND_CH / 3.0)          # cos(109.5) = -1/3
            r = BOND_CH * (2 * np.sqrt(2) / 3.0)        # sin(109.5)
            hydro.append(Atom(f'H{i+1}A', base + p2 * r, 'hc', 0.06))
            hydro.append(Atom(f'H{i+1}B', base + (-0.5 * p2 + 0.866 * p1) * r, 'hc', 0.06))
            hydro.append(Atom(f'H{i+1}C', base + (-0.5 * p2 - 0.866 * p1) * r, 'hc', 0.06))
            continue
        if i == 0:
            v1 = heavy[0].pos - cp          # toward distal carboxyl C
            v2 = c_pos[1] - cp
        elif i == n_chain - 1:              # diacid only: neighbor is C1P
            v1 = c_pos[i-1] - cp
            v2 = heavy[3 + n_chain].pos - cp
        else:
            v1 = c_pos[i-1] - cp
            v2 = c_pos[i+1] - cp
        add_h(hydro, f'H{i+1}', cp, v1, v2)

    all_atoms = heavy + hydro
    target_q = -1.0 if kind == 'mono' else -2.0
    total_q = sum(a.charge for a in all_atoms)
    for a in all_atoms:
        a.charge *= target_q / total_q
    # Round to mol2's 4 decimals; residual onto a mid-chain carbon so the
    # written charges sum EXACTLY to the integer target.
    for a in all_atoms:
        a.charge = round(a.charge, 4)
    resid = target_q - sum(a.charge for a in all_atoms)
    anchor_name = 'C09' if n_chain >= 9 else f'C{n_chain:02d}'
    mid_c = next(a for a in all_atoms if a.name == anchor_name)
    mid_c.charge = round(mid_c.charge + resid, 4)

    # ---- Bonds ----
    bonds = []
    # heavy indices (1-based): 1=C1X/C1D, 2=O1D, 3=O2D, 4..3+n_chain=C01..Cn
    bonds.append((1, 2, 2)); bonds.append((1, 3, 1)); bonds.append((1, 4, 1))
    for i in range(n_chain - 1):
        bonds.append((4 + i, 5 + i, 1))
    if kind == 'diacid':
        c1p = 4 + n_chain          # 1-based index of C1P
        bonds.append((c1p - 1, c1p, 1)); bonds.append((c1p, c1p + 1, 2)); bonds.append((c1p, c1p + 2, 1))
    n_heavy = len(heavy)
    for h_idx, a in enumerate(all_atoms[n_heavy:]):
        cnum = int(a.name[1:].rstrip('ABC'))
        bonds.append((4 + cnum - 1, n_heavy + h_idx + 1, 1))

    # ---- Write mol2 ----
    out = os.path.join(HERE, f"{label}_fa3.mol2")
    with open(out, 'w') as f:
        f.write("@<TRIPOS>MOLECULE\nFAH\n")
        f.write(f" {len(all_atoms)} {len(bonds)} 1 0 0\nSMALL\nGAFF2\n"
                f"{label} at FA3 (all-atom)\n\n")
        f.write("@<TRIPOS>ATOM\n")
        for aid, a in enumerate(all_atoms):
            f.write(f"{aid+1:6d} {a.name:6s} {a.pos[0]:10.4f} {a.pos[1]:10.4f} {a.pos[2]:10.4f} "
                    f"{a.atype:4s} 1 FAH {a.charge:10.4f}\n")
        f.write("@<TRIPOS>BOND\n")
        for bid, (i, j, t) in enumerate(bonds):
            f.write(f"{bid+1:6d} {i:6d} {j:6d} {t}\n")
        f.write("@<TRIPOS>SUBSTRUCTURE\n1 FAH 1 TEMP 0 **** **** 0 ROOT\n")

    # ---- Report ----
    o1d, o2d = heavy[1].pos, heavy[2].pos
    print(f"[{label}] wrote {out}")
    print(f"  atoms={len(all_atoms)} (heavy={n_heavy}, H={len(hydro)}), bonds={len(bonds)}, "
          f"charge={sum(a.charge for a in all_atoms):+.3f}")
    print(f"  O1D-O2D: {np.linalg.norm(o1d-o2d):.2f} A | "
          f"O1D->MYR O1: {np.linalg.norm(o1d-MYR['O1']):.2f} A, O2D->MYR O2: {np.linalg.norm(o2d-MYR['O2']):.2f} A")
    # self-overlap guard: no two heavy atoms closer than 1.2 A unless bonded
    bonded = {tuple(sorted((i, j))) for i, j, _ in bonds}
    worst = min((np.linalg.norm(all_atoms[i-1].pos - all_atoms[j-1].pos), i, j)
                for i in range(1, n_heavy + 1) for j in range(i + 1, n_heavy + 1)
                if (i, j) not in bonded)
    flag = "OK" if worst[0] > 1.2 else "FAIL -- CHAIN SELF-OVERLAP"
    print(f"  min nonbonded heavy-heavy distance: {worst[0]:.2f} A [{flag}]")
    # clash guard vs HSA
    fah = np.array([a.pos for a in heavy])
    dmin = np.linalg.norm(fah[:, None, :] - HSA[None, :, :], axis=-1).min()
    print(f"  min FA-heavy..HSA-heavy distance: {dmin:.2f} A (want > ~1.5)")
    if kind == 'diacid':
        o1p, o2p = heavy[-2].pos, heavy[-1].pos
        print(f"  O1P-O2P: {np.linalg.norm(o1p-o2p):.2f} A (carboxyl, ~2.23)")

if __name__ == "__main__":
    args = sys.argv[1:]
    if not args or args[0] not in ('mono', 'diacid', 'both'):
        print(__doc__)
        sys.exit(1)
    kind = args[0]
    force = '--force' in args
    label_override = None
    if '--label' in args:
        i = args.index('--label')
        if i + 1 >= len(args) or args[i + 1].startswith('-'):
            sys.exit("--label needs a value")
        label_override = args[i + 1]
    lengths = [int(x) for x in args[1:] if x.lstrip('-').isdigit()]
    if not lengths:
        print(__doc__)
        sys.exit(1)
    kinds = ['mono', 'diacid'] if kind == 'both' else [kind]
    if label_override and len(lengths) * len(kinds) != 1:
        sys.exit("--label only valid with exactly one system")
    for L in lengths:
        for k in kinds:
            n_chain = L - 2 if k == 'diacid' else L - 1
            label = label_override or f"c{L}_{'diacid' if k == 'diacid' else 'monoacid'}"
            out = os.path.join(HERE, f"{label}_fa3.mol2")
            if os.path.exists(out) and not force:
                print(f"SKIP {label}: {out} exists (use --force to overwrite)")
                continue
            build(k, n_chain, label)
