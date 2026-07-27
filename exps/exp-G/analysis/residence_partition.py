#!/usr/bin/env python3
"""exp-G residence-partition analysis (2026-07-25).

Reads tleap/complex.prmtop + md/rep{1,2,3}/traj.dcd (1003 frames x 100 ps,
dt derived from the DCDReporter interval of 50,000 steps x 2 fs -- file
timestamps are not trusted). The three replicas (seeds 101/202/303) are
independent trajectories: all statistics are computed per rep first, then
merged as mean(rep-means) +/- SD over n=3 reps.

Per-frame metrics (minimum-image distances, orthorhombic box per frame):
  d(O56/O57 -> ARG346/ARG483 guanidinium), tail(O56/57)->HSA nearest heavy,
  d(O38 -> ARG408 guanidinium), peptide-ECD min heavy distance,
  ECD-HSA COM distance (CA-based), LNK end-to-end (C11->C55), LNK Rg.

Residence classes for the distal carboxyl:
  FA3      : min(dO56_R346, dO56_R483, dO57_R346, dO57_R483) < 4 A
  surface  : not FA3 and tail->HSA < 4 A
  aqueous  : tail->HSA > 10 A
  transition: 4-10 A (listed separately)

Events:
  de-anchor : all 4 anchor distances > 6 A sustained for >= 10 frames (>=1 ns)
  detach    : peptide-ECD min distance > 6 A sustained >= 10 frames

Outputs: frames_rep{N}.csv, summary.json (used to write RESULTS.md).
Also performs the rep-1 frame-0 bond-length check (trajectory-topology match).
"""
import os, json
import numpy as np
import parmed as pmd
import MDAnalysis as mda
import warnings
warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))
EXP_G = os.path.normpath(os.path.join(HERE, ".."))
REPS = [1, 2, 3]
DT_PS = 100.0
SUSTAIN = 10  # frames = 1.0 ns

# ---------------- topology groups (1-based residue positions) ----------------
parm = pmd.load_file(os.path.join(EXP_G, "tleap", "complex.prmtop"))
prot_res = [r for r in parm.residues if r.name not in ("WAT", "HOH", "Na+", "Cl-", "NA", "CL")]
pos_of = {id(r): i + 1 for i, r in enumerate(prot_res)}
assert len(prot_res) == 709

def rpos(a):
    return pos_of.get(id(a.residue), 10 ** 9)  # waters/ions -> +inf

hsa_heavy = [a.idx for a in parm.atoms if rpos(a) <= 582 and a.element != 1]
ecd_heavy = [a.idx for a in parm.atoms if 583 <= rpos(a) <= 682 and a.element != 1]
pep_heavy = [a.idx for a in parm.atoms if 683 <= rpos(a) <= 708 and a.element != 1]
hsa_ca = [a.idx for a in parm.atoms if rpos(a) <= 582 and a.name == "CA"]
ecd_ca = [a.idx for a in parm.atoms if 583 <= rpos(a) <= 682 and a.name == "CA"]
lnk = {a.name: a.idx for a in parm.atoms if a.residue.name == "LNK"}
lnk_heavy = [a.idx for a in parm.atoms if a.residue.name == "LNK" and a.element != 1]

def guanid(rnum):
    r = prot_res[rnum - 1]
    assert r.name == "ARG", (rnum, r.name)
    return [a.idx for a in r.atoms if a.name in ("NE", "NH1", "NH2")]

R346, R483, R408 = guanid(346), guanid(483), guanid(408)
nz_idx = next(a.idx for a in parm.atoms
              if a.name == "NZ" and a.residue.name == "LYS" and rpos(a) == 699)
print(f"groups: HSA {len(hsa_heavy)}, ECD {len(ecd_heavy)}, pep {len(pep_heavy)}, "
      f"LNK {len(lnk_heavy)}, ARG346/483/408 guanidinium ok")

def mdist(A, B, box):
    d = A[:, None, :] - B[None, :, :]
    d -= box[:3] * np.round(d / box[:3])
    return float(np.sqrt((d ** 2).sum(-1)).min())

def rg(x):
    c = x.mean(axis=0)
    return float(np.sqrt(((x - c) ** 2).sum(-1).mean()))

# ---------------- frame-0 bond check (rep1) ----------------
u = mda.Universe(os.path.join(EXP_G, "tleap", "complex.prmtop"),
                 os.path.join(EXP_G, "md", "rep1", "traj.dcd"))
ts = u.trajectory[0]
xyz, box = ts.positions, ts.dimensions
checks = {
    "NZ-C11 (1.38)": np.linalg.norm(xyz[nz_idx] - xyz[lnk["C11"]]),
    "C55-O56 (1.26)": np.linalg.norm(xyz[lnk["C55"]] - xyz[lnk["O56"]]),
    "C37-O38 (1.24)": np.linalg.norm(xyz[lnk["C37"]] - xyz[lnk["O38"]]),
    "O56->ARG483 (2.78)": mdist(xyz[[lnk["O56"]]], xyz[R483], box),
    "O38->ARG408 (3.95)": mdist(xyz[[lnk["O38"]]], xyz[R408], box),
}
print("rep1 frame-0 bond/anchor check:")
for k, v in checks.items():
    print(f"  {k}: {v:.2f} A")
assert abs(checks["NZ-C11 (1.38)"] - 1.38) < 0.2, "trajectory-topology mismatch!"
print("  -> trajectory matches topology\n")

# ---------------- per-frame metric loop ----------------
def analyze(rep):
    u = mda.Universe(os.path.join(EXP_G, "tleap", "complex.prmtop"),
                     os.path.join(EXP_G, "md", f"rep{rep}", "traj.dcd"))
    rows = []
    for ts in u.trajectory:
        xyz, box = ts.positions, ts.dimensions
        o56, o57 = xyz[[lnk["O56"]]], xyz[[lnk["O57"]]]
        d56_346 = mdist(o56, xyz[R346], box); d56_483 = mdist(o56, xyz[R483], box)
        d57_346 = mdist(o57, xyz[R346], box); d57_483 = mdist(o57, xyz[R483], box)
        anchor = min(d56_346, d56_483, d57_346, d57_483)
        tail = min(mdist(o56, xyz[hsa_heavy], box), mdist(o57, xyz[hsa_heavy], box))
        o38 = mdist(xyz[[lnk["O38"]]], xyz[R408], box)
        if anchor < 4.0: cls = "FA3"
        elif tail < 4.0: cls = "surface"
        elif tail > 10.0: cls = "aqueous"
        else: cls = "transition"
        pep_ecd = mdist(xyz[pep_heavy], xyz[ecd_heavy], box)
        com = float(np.linalg.norm(xyz[hsa_ca].mean(axis=0) - xyz[ecd_ca].mean(axis=0)))
        e2e = float(np.linalg.norm(xyz[lnk["C11"]] - xyz[lnk["C55"]]))
        rows.append(dict(frame=ts.frame, time_ps=ts.frame * DT_PS,
                         dO56_R346=d56_346, dO56_R483=d56_483,
                         dO57_R346=d57_346, dO57_R483=d57_483,
                         anchor_min=anchor, tail_hsa_min=tail, cls=cls,
                         dO38_R408=o38, pep_ecd_min=pep_ecd,
                         com_hsa_ecd=com, lnk_e2e=e2e, lnk_rg=rg(xyz[lnk_heavy])))
    return rows

def first_sustained(cond, sustain=SUSTAIN):
    """First index where cond holds for all of the next `sustain` frames."""
    c = np.asarray(cond, dtype=bool)
    run = np.concatenate(([0], np.cumsum(c)))
    for i in range(len(c) - sustain + 1):
        if run[i + sustain] - run[i] == sustain:
            return i
    return None

summary = {}
all_rows = {}
for rep in REPS:
    rows = analyze(rep)
    all_rows[rep] = rows
    # CSV
    csv_path = os.path.join(HERE, f"frames_rep{rep}.csv")
    with open(csv_path, "w") as f:
        f.write(",".join(rows[0].keys()) + "\n")
        for r in rows:
            f.write(",".join(f"{r[k]:.3f}" if isinstance(r[k], float) else str(r[k]) for k in r) + "\n")
    n = len(rows)
    fracs = {c: sum(1 for r in rows if r["cls"] == c) / n for c in ("FA3", "surface", "transition", "aqueous")}
    anchor_all6 = [r["anchor_min"] > 6.0 for r in rows]
    # de-anchor requires BOTH O's away from BOTH Args: anchor_min>6 covers it
    i_de = first_sustained(anchor_all6)
    i_dt = first_sustained([r["pep_ecd_min"] > 6.0 for r in rows])
    tail = np.array([r["tail_hsa_min"] for r in rows])
    pep = np.array([r["pep_ecd_min"] for r in rows])
    summary[rep] = dict(
        frames=n, ns=n * DT_PS / 1000,
        frac=fracs,
        de_anchor_ps=(rows[i_de]["time_ps"] if i_de is not None else None),
        reanchor_after_de=bool(i_de is not None and any(r["cls"] == "FA3" for r in rows[i_de + SUSTAIN:])),
        o38_bridge_frac=sum(1 for r in rows if r["dO38_R408"] < 4.0) / n,
        o38_mean=float(np.mean([r["dO38_R408"] for r in rows])),
        pep_ecd_mean=float(pep.mean()), pep_ecd_sd=float(pep.std()),
        detach_ps=(rows[i_dt]["time_ps"] if i_dt is not None else None),
        tail_mean=float(tail.mean()), tail_sd=float(tail.std()),
        tail_last25_mean=float(tail[-n // 4:].mean()), tail_last25_sd=float(tail[-n // 4:].std()),
        pep_last25_mean=float(pep[-n // 4:].mean()), pep_last25_sd=float(pep[-n // 4:].std()),
        com_mean=float(np.mean([r["com_hsa_ecd"] for r in rows])),
        com_sd=float(np.std([r["com_hsa_ecd"] for r in rows])),
        e2e_mean=float(np.mean([r["lnk_e2e"] for r in rows])),
        e2e_sd=float(np.std([r["lnk_e2e"] for r in rows])),
        rg_mean=float(np.mean([r["lnk_rg"] for r in rows])),
        rg_sd=float(np.std([r["lnk_rg"] for r in rows])),
    )
    print(f"rep{rep}: {n} frames, FA3 {fracs['FA3']:.1%}, surface {fracs['surface']:.1%}, "
          f"transition {fracs['transition']:.1%}, aqueous {fracs['aqueous']:.1%}, "
          f"de-anchor@{summary[rep]['de_anchor_ps']} ps, O38 bridge {summary[rep]['o38_bridge_frac']:.1%}")

# merged: mean +/- SD over rep-means (n=3)
merged = {}
for key in ("o38_bridge_frac", "pep_ecd_mean", "pep_ecd_sd", "com_mean", "com_sd",
            "e2e_mean", "e2e_sd", "rg_mean", "rg_sd", "o38_mean",
            "tail_mean", "tail_sd", "tail_last25_mean", "tail_last25_sd",
            "pep_last25_mean", "pep_last25_sd"):
    vals = [summary[r][key] for r in REPS]
    merged[key] = (float(np.mean(vals)), float(np.std(vals)))
for c in ("FA3", "surface", "transition", "aqueous"):
    vals = [summary[r]["frac"][c] for r in REPS]
    merged.setdefault("frac", {})[c] = (float(np.mean(vals)), float(np.std(vals)))
summary["merged"] = merged

with open(os.path.join(HERE, "summary.json"), "w") as f:
    json.dump(summary, f, indent=2, default=str)
print(f"\nwrote frames_rep{{1,2,3}}.csv + summary.json")
