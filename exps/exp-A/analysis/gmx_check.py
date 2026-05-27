#!/usr/bin/env python3
"""Analyze GROMACS trajectories: RMSD, catalytic distance, WT vs Aib8 comparison."""
import mdtraj as md
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

REPO = "/home/scroll/personal/semaglutide-in-silico"
OUT = f"{REPO}/exps/exp-A/analysis"
os.makedirs(OUT, exist_ok=True)

for label in ["wt", "aib8"]:
    gro = f"{REPO}/exps/exp-A/gmx/{label}/npt.gro"
    xtc = f"{REPO}/exps/exp-A/gmx/{label}/md.xtc"
    if not os.path.exists(xtc):
        print(f"{label}: no xtc"); continue

    print(f"\n=== GROMACS {label.upper()} ===")
    t = md.load(xtc, top=gro)
    dt = 0.1  # 100 ps per frame (nstxout-compressed=50000 * 2fs)
    time_ns = np.arange(t.n_frames) * dt
    half = t.n_frames // 2
    print(f"  {t.n_frames} frames, {t.n_atoms} atoms, {time_ns[-1]:.1f} ns")

    # In GROMACS gro (from ParmEd): residues are numbered sequentially
    # DPP-4 is residues 1-728, peptide is 729-759
    # Find them by looking for the peptide start (HIS/HIE near the end)
    residues = list(t.topology.residues)
    pep_start = None
    for r in reversed(residues):
        if r.name in ('HIE', 'HIS') and r.resSeq > 700:
            pep_start = r.resSeq
            break

    if pep_start is None:
        # Fallback: last non-water residue minus ~31
        non_wat = [r for r in residues if r.name not in ('SOL','WAT','NA','CL','Na+','Cl-')]
        pep_start = non_wat[-1].resSeq - 30

    pep_end = pep_start + 30
    dpp4_start = 1
    dpp4_end = pep_start - 1
    print(f"  DPP-4: resid {dpp4_start}–{dpp4_end}, Peptide: resid {pep_start}–{pep_end}")

    # Select atoms
    dpp4_bb = t.topology.select(f"resid {dpp4_start} to {dpp4_end} and backbone")
    pep_bb = t.topology.select(f"resid {pep_start} to {pep_end} and backbone")

    # Find catalytic Ser OG (closest to peptide N-term)
    all_ser_og = t.topology.select("resname SER and name OG")
    pep_nterm_ca = t.topology.select(f"resid {pep_start} to {pep_start+3} and name CA")
    if len(all_ser_og) > 0 and len(pep_nterm_ca) > 0:
        d = np.linalg.norm(t.xyz[0, all_ser_og] - t.xyz[0, pep_nterm_ca].mean(0), axis=1)
        ser_og_idx = all_ser_og[np.argmin(d)]
        ser_res = t.topology.atom(ser_og_idx).residue
        print(f"  Catalytic Ser: {ser_res.name}{ser_res.resSeq} OG (dist={d.min():.1f} Å)")

    # Peptide residue 2 C (Ala8 or Aib8)
    pep2_c = t.topology.select(f"resid {pep_start+1} and name C")
    if len(pep2_c) == 0:
        pep2_c = t.topology.select(f"(resname ALA or resname AIB) and name C and resid > {pep_start-5}")
    if len(pep2_c) > 0:
        r = t.topology.atom(pep2_c[0]).residue
        print(f"  Peptide res2: {r.name}{r.resSeq} C")

    if len(dpp4_bb) == 0 or len(pep_bb) == 0:
        print("  SKIP: empty selection"); continue

    # Align on DPP-4 and compute
    t.superpose(t[0], atom_indices=dpp4_bb)
    dpp4_rmsd = md.rmsd(t, t[0], atom_indices=dpp4_bb) * 10
    pep_rmsd = md.rmsd(t, t[0], atom_indices=pep_bb) * 10
    cat_dist = np.linalg.norm(t.xyz[:, ser_og_idx] - t.xyz[:, pep2_c[0]], axis=1) * 10

    print(f"  DPP-4 RMSD (2nd half): {dpp4_rmsd[half:].mean():.1f} ± {dpp4_rmsd[half:].std():.1f} Å")
    print(f"  Peptide RMSD (2nd half): {pep_rmsd[half:].mean():.1f} ± {pep_rmsd[half:].std():.1f} Å")
    print(f"  Catalytic dist (2nd half): {cat_dist[half:].mean():.1f} ± {cat_dist[half:].std():.1f} Å")
    print(f"  Dist range: {cat_dist.min():.1f}–{cat_dist.max():.1f} Å")

    # Find equilibrium (last 50 ns)
    last_50 = int(50 / dt)  # frames in last 50 ns
    if last_50 > 0 and t.n_frames > last_50:
        eq_start = t.n_frames - last_50
        print(f"  Last 50 ns: DPP4={dpp4_rmsd[eq_start:].mean():.1f}±{dpp4_rmsd[eq_start:].std():.1f}Å  "
              f"Pep={pep_rmsd[eq_start:].mean():.1f}±{pep_rmsd[eq_start:].std():.1f}Å  "
              f"CatD={cat_dist[eq_start:].mean():.1f}±{cat_dist[eq_start:].std():.1f}Å")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes[0,0].plot(time_ns, dpp4_rmsd, alpha=0.7); axes[0,0].set_ylabel("Å")
    axes[0,0].set_title(f"DPP-4 backbone RMSD ({dpp4_rmsd.mean():.1f} Å)")
    axes[0,1].plot(time_ns, pep_rmsd, alpha=0.7); axes[0,1].set_ylabel("Å")
    axes[0,1].set_title(f"Peptide backbone RMSD ({pep_rmsd.mean():.1f} Å)")
    axes[1,0].plot(time_ns, cat_dist, alpha=0.7, color='red')
    axes[1,0].axhline(y=3.5, color='gray', linestyle='--', label='ref ~3.5Å')
    axes[1,0].set_ylabel("Å"); axes[1,0].set_title(f"Ser OG → Peptide Res2 C"); axes[1,0].legend()
    axes[1,1].hist(cat_dist[half:], bins=50, alpha=0.7)
    axes[1,1].axvline(x=3.5, color='gray', linestyle='--')
    axes[1,1].set_title("Catalytic dist distribution (2nd half)")
    for ax in axes.flat: ax.set_xlabel("Time (ns)")
    plt.suptitle(f"GROMACS {label.upper()}: {time_ns[-1]:.0f} ns", fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{OUT}/gmx_{label}_quick.png", dpi=100)
    plt.close()
    print(f"  Plot: {OUT}/gmx_{label}_quick.png")

print("\nDone!")
