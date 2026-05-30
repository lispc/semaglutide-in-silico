#!/usr/bin/env python3
"""Deep analysis of ACE-linker MD: FA3 binding stability, anchoring, tail dynamics."""
import mdtraj as md, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os, sys

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_C = f"{REPO}/exps/exp-C"
prmtop = f"{EXP_C}/tleap/linker_ace.prmtop"
md_dir = f"{EXP_C}/md/linker_ace"
OUT = f"{EXP_C}/analysis"

# --- Load all 3 replicas ---
trajs = {}
for rep in [1, 2, 3]:
    dcd = f"{md_dir}/rep{rep}/linker_ace_traj.dcd"
    if not os.path.exists(dcd):
        print(f"Missing: {dcd}")
        continue
    # Use iterload to handle incomplete DCD (MD may still be writing)
    frames = []
    for chunk in md.iterload(dcd, top=prmtop, chunk=100):
        frames.append(chunk)
        if chunk.n_frames < 100:
            break  # last chunk
    if not frames:
        print(f"Rep {rep}: no frames read, skipping")
        continue
    t = frames[0] if len(frames) == 1 else md.join(frames)
    print(f"Rep {rep}: {t.n_frames} frames, {t.time[-1]:.0f} ps")
    trajs[rep] = t

if not trajs:
    print("No trajectories found!")
    sys.exit(1)

# --- Atom selections ---
t0 = list(trajs.values())[0]
hsa_ca = t0.topology.select("protein and name CA")
print(f"HSA CA: {len(hsa_ca)}")

# LFA heavy atoms
lfa_heavy = t0.topology.select("resname LFA and not element H")
print(f"LFA heavy: {len(lfa_heavy)}")

# C18 tail: CH2 chain + terminal carboxylate (C30-C46 + O47,O48)
# The tail carbons are the ones after the last OEG amide through the terminal COO-
# From mol2: C30(idx 9249) through C46(idx 9265), O47(idx 9266), O48(idx 9267)
tail_c = t0.topology.select("resname LFA and (name C30 C31 C32 C33 C34 C35 C36 C37 C38 C39 C40 C41 C42 C43 C44 C45 C46)")
tail_heavy = t0.topology.select("resname LFA and (name C30 C31 C32 C33 C34 C35 C36 C37 C38 C39 C40 C41 C42 C43 C44 C45 C46 O47 O48)")
distal_o = t0.topology.select("resname LFA and (name O47 O48)")
proximal_o = t0.topology.select("resname LFA and (name O50 O51)")

# ARG482 guanidinium
arg482_n = t0.topology.select("(resSeq 482 or resid 482) and (name NH1 NH2)")
arg348_n = t0.topology.select("(resSeq 348 or resid 348) and (name NH1 NH2)")
print(f"Tail heavy: {len(tail_heavy)}, Distal O: {len(distal_o)}, ARG482 N: {len(arg482_n)}, ARG348 N: {len(arg348_n)}")

# If ARG482 selection failed by resSeq, try by index
if len(arg482_n) == 0:
    # ARG482 NH1=7668, NH2=7671
    arg482_n = np.array([7668, 7671])
    print(f"  ARG482 fallback: {arg482_n}")

# --- Analysis ---
results = {}
for rep, t in trajs.items():
    # Align to first frame
    t.superpose(t[0], atom_indices=hsa_ca)

    # HSA backbone RMSD
    hsa_rmsd = md.rmsd(t, t[0], atom_indices=hsa_ca) * 10  # nm→Å

    # C18 tail RMSD
    tail_rmsd = md.rmsd(t, t[0], atom_indices=tail_heavy) * 10

    # Distal COO⁻ → ARG482 (minimum distance between any O..any N)
    d_482 = []
    for f in range(t.n_frames):
        od_xyz = t.xyz[f, distal_o]
        ar_xyz = t.xyz[f, arg482_n]
        dmin = min(np.linalg.norm(o - n) for o in od_xyz for n in ar_xyz)
        d_482.append(dmin * 10)  # nm→Å
    d_482 = np.array(d_482)

    # Distal COO⁻ → ARG348 (second anchor candidate)
    d_348 = []
    if len(arg348_n) > 0:
        for f in range(t.n_frames):
            od_xyz = t.xyz[f, distal_o]
            ar_xyz = t.xyz[f, arg348_n]
            dmin = min(np.linalg.norm(o - n) for o in od_xyz for n in ar_xyz)
            d_348.append(dmin * 10)
        d_348 = np.array(d_348)
    else:
        d_348 = np.zeros(t.n_frames)

    # FA COM displacement (C18 tail COM vs frame 0)
    tail_com_0 = t.xyz[0, tail_c].mean(axis=0)
    tail_com = t.xyz[:, tail_c].mean(axis=1)
    com_disp = np.linalg.norm(tail_com - tail_com_0, axis=1) * 10  # nm→Å

    # Proximal COO⁻ → nearest basic residue
    d_prox_482 = []
    if len(proximal_o) > 0:
        for f in range(t.n_frames):
            op_xyz = t.xyz[f, proximal_o]
            ar_xyz = t.xyz[f, arg482_n]
            dmin = min(np.linalg.norm(o - n) for o in op_xyz for n in ar_xyz)
            d_prox_482.append(dmin * 10)
        d_prox_482 = np.array(d_prox_482)
    else:
        d_prox_482 = np.zeros(t.n_frames)

    # End-to-end distance of linker (ACE C → distal COO⁻ C, C1→C46)
    c1_idx = t0.topology.select("resname LFA and name C1")
    c46_idx = t0.topology.select("resname LFA and name C46")
    if len(c1_idx) > 0 and len(c46_idx) > 0:
        ee = np.linalg.norm(t.xyz[:, c46_idx[0]] - t.xyz[:, c1_idx[0]], axis=1) * 10
    else:
        ee = np.zeros(t.n_frames)

    # Halfway point for "equilibrated" stats
    half = t.n_frames // 2

    results[rep] = dict(hsa_rmsd=hsa_rmsd, tail_rmsd=tail_rmsd, d_482=d_482,
                        d_348=d_348, com_disp=com_disp, d_prox_482=d_prox_482,
                        ee=ee, half=half, n_frames=t.n_frames)

    print(f"\nRep {rep} (2nd half avg):")
    print(f"  HSA RMSD:        {hsa_rmsd[half:].mean():.2f}±{hsa_rmsd[half:].std():.2f} Å")
    print(f"  Tail RMSD:       {tail_rmsd[half:].mean():.2f}±{tail_rmsd[half:].std():.2f} Å")
    print(f"  COO⁻→ARG482:     {d_482[half:].mean():.2f}±{d_482[half:].std():.2f} Å")
    print(f"  COO⁻→ARG348:     {d_348[half:].mean():.2f}±{d_348[half:].std():.2f} Å")
    print(f"  Tail COM displ:  {com_disp[half:].mean():.2f}±{com_disp[half:].std():.2f} Å")
    print(f"  ProxCOO⁻→ARG482: {d_prox_482[half:].mean():.2f}±{d_prox_482[half:].std():.2f} Å")
    print(f"  Linker E2E:      {ee[half:].mean():.2f}±{ee[half:].std():.2f} Å")

    # Also: what fraction of frames is d_482 < 5 Å (bound)?
    frac_bound = (d_482[half:] < 5.0).mean() * 100
    print(f"  Fraction COO⁻-ARG482 <5Å: {frac_bound:.1f}%")
    # Fraction d_482 < 10 Å (proximal)
    frac_near = (d_482[half:] < 10.0).mean() * 100
    print(f"  Fraction COO⁻-ARG482 <10Å: {frac_near:.1f}%")

# === Time series plot ===
fig, axes = plt.subplots(2, 3, figsize=(20, 10))
rep_colors = {1: '#2196F3', 2: '#4CAF50', 3: '#FF9800'}
time_labels = {1: [], 2: [], 3: []}  # store time arrays per rep

for rep in [1, 2, 3]:
    if rep not in results:
        continue
    r = results[rep]
    time_ns = np.arange(r['n_frames']) * (50000 * 0.002 * 1e-3)  # 50k steps/frame * 2fs * 1e-3ps
    c = rep_colors[rep]

    axes[0, 0].plot(time_ns, r['hsa_rmsd'], color=c, alpha=0.7, label=f'Rep {rep}')
    axes[0, 1].plot(time_ns, r['tail_rmsd'], color=c, alpha=0.7, label=f'Rep {rep}')
    axes[0, 2].plot(time_ns, r['d_482'], color=c, alpha=0.7, label=f'Rep {rep}')
    axes[1, 0].plot(time_ns, r['com_disp'], color=c, alpha=0.7, label=f'Rep {rep}')
    axes[1, 1].plot(time_ns, r['ee'], color=c, alpha=0.7, label=f'Rep {rep}')
    axes[1, 2].plot(time_ns, r['d_348'], color=c, alpha=0.3, label=f'Rep {rep}')

axes[0, 0].set_ylabel('HSA CA RMSD (Å)')
axes[0, 0].set_title('HSA Backbone Stability')
axes[0, 1].set_ylabel('C18 Tail RMSD (Å)')
axes[0, 1].set_title('C18 Fatty Acid Tail RMSD')
axes[0, 2].set_ylabel('Distance (Å)')
axes[0, 2].set_title('Distal COO⁻ → ARG482')
axes[0, 2].axhline(3.5, color='gray', linestyle='--', alpha=0.5, label='H-bond cutoff')
axes[1, 0].set_ylabel('Displacement (Å)')
axes[1, 0].set_title('C18 Tail COM Displacement')
axes[1, 1].set_ylabel('Distance (Å)')
axes[1, 1].set_title('Linker End-to-End (ACE-C→C18-COO⁻)')
axes[1, 2].set_ylabel('Distance (Å)')
axes[1, 2].set_title('Distal COO⁻ → ARG348')

for ax in axes.flat:
    ax.legend(fontsize=7, loc='upper right')
    ax.set_xlabel('Time (ns)')

plt.tight_layout()
plt.savefig(f"{OUT}/ace_linker_analysis.png", dpi=150)
print(f"\nSaved: {OUT}/ace_linker_analysis.png")

# === Summary table ===
print(f"\n{'='*90}")
print(f"  {'Rep':<5s} {'HSA RMSD':>12s} {'Tail RMSD':>12s} {'482 dist':>10s} {'348 dist':>10s} {'COM disp':>10s} {'E2E':>10s} {'<5Å frac':>10s}")
print(f"  {'-'*88}")
for rep in [1, 2, 3]:
    if rep not in results:
        continue
    r = results[rep]
    h = r['half']
    hsa = f"{r['hsa_rmsd'][h:].mean():.1f}±{r['hsa_rmsd'][h:].std():.1f}"
    tail = f"{r['tail_rmsd'][h:].mean():.1f}±{r['tail_rmsd'][h:].std():.1f}"
    d482 = f"{r['d_482'][h:].mean():.1f}±{r['d_482'][h:].std():.1f}"
    d348 = f"{r['d_348'][h:].mean():.1f}±{r['d_348'][h:].std():.1f}"
    com = f"{r['com_disp'][h:].mean():.1f}±{r['com_disp'][h:].std():.1f}"
    e2e = f"{r['ee'][h:].mean():.1f}±{r['ee'][h:].std():.1f}"
    frac = f"{(r['d_482'][h:] < 5.0).mean()*100:.1f}%"
    print(f"  {rep:<5d} {hsa:>12s} {tail:>12s} {d482:>10s} {d348:>10s} {com:>10s} {e2e:>10s} {frac:>10s}")
