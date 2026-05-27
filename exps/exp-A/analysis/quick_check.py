#!/usr/bin/env python3
"""Quick MD analysis: RMSD, catalytic distance, basic health check."""
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
    if label == "wt":
        prmtop = f"{REPO}/exps/exp-A/tleap/wt.prmtop"
    else:
        prmtop = f"{REPO}/exps/exp-A/md/aib8/aib8_modified.prmtop"
    dcd = f"{REPO}/exps/exp-A/md/{label}/{label}_traj.dcd"
    if not os.path.exists(dcd): print(f"{label}: no DCD"); continue

    print(f"\n=== {label.upper()} ===")
    t = md.load(dcd, top=prmtop)
    dt_ns = 0.1  # 100 ps per frame
    time_ns = np.arange(t.n_frames) * dt_ns
    half = t.n_frames // 2
    print(f"  {t.n_frames} frames, {time_ns[-1]:.1f} ns")

    # DPP-4 backbone (resid 0-727) for RMSD alignment
    dpp4_idx = t.topology.select("resid 0 to 727 and backbone")
    # Peptide backbone (resid 728-758)
    pep_idx = t.topology.select("resid 728 to 758 and backbone")
    # Ser630 OG
    ser_og = t.topology.select("resid 630 and name OG")
    if len(ser_og) == 0:
        # tleap may have shifted numbering; find closest SER OG to peptide
        all_ser = t.topology.select("resname SER and name OG")
        pep_near = t.topology.select("resid 728 to 733 and name CA")
        if len(all_ser) > 0 and len(pep_near) > 0:
            d = np.linalg.norm(t.xyz[0,all_ser]-t.xyz[0,pep_near].mean(0), axis=1)
            ser_og = [all_ser[np.argmin(d)]]
            r = t.topology.atom(ser_og[0]).residue
            print(f"  Ser OG: resid {r.resSeq} {r.name} (auto-detected)")
    # Ala8/Aib8 C (resid 729)
    pep2_c = t.topology.select("resid 729 and name C")
    print(f"  Pep res2 C: {len(pep2_c)}, Ser OG: {len(ser_og)}")

    if len(ser_og) == 0 or len(pep2_c) == 0: print("  SKIP"); continue

    # Align all frames to frame 0 on DPP-4 backbone, compute RMSD
    t.superpose(t[0], atom_indices=dpp4_idx)
    dpp4_rmsd = md.rmsd(t, t[0], atom_indices=dpp4_idx) * 10  # nm→Å
    pep_rmsd = md.rmsd(t, t[0], atom_indices=pep_idx) * 10 if len(pep_idx) > 0 else np.zeros(t.n_frames)

    # Catalytic distance
    cat_dist = np.linalg.norm(t.xyz[:, ser_og[0]] - t.xyz[:, pep2_c[0]], axis=1) * 10

    print(f"  DPP-4 RMSD (2nd half): {dpp4_rmsd[half:].mean():.1f} ± {dpp4_rmsd[half:].std():.1f} Å")
    print(f"  Peptide RMSD (2nd half): {pep_rmsd[half:].mean():.1f} ± {pep_rmsd[half:].std():.1f} Å")
    print(f"  Catalytic dist (2nd half): {cat_dist[half:].mean():.1f} ± {cat_dist[half:].std():.1f} Å")
    print(f"  Catalytic dist range: {cat_dist.min():.1f}–{cat_dist.max():.1f} Å")

    # Plot
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes[0,0].plot(time_ns, dpp4_rmsd, alpha=0.7); axes[0,0].set_ylabel("Å")
    axes[0,0].set_title(f"DPP-4 backbone RMSD (mean={dpp4_rmsd.mean():.1f}Å)")
    axes[0,1].plot(time_ns, pep_rmsd, alpha=0.7); axes[0,1].set_ylabel("Å")
    axes[0,1].set_title(f"Peptide backbone RMSD (mean={pep_rmsd.mean():.1f}Å)")
    axes[1,0].plot(time_ns, cat_dist, alpha=0.7, color='red')
    axes[1,0].axhline(y=3.5, color='gray', linestyle='--', label='ref ~3.5Å')
    axes[1,0].set_ylabel("Å"); axes[1,0].set_title("Ser630 OG → Ala8/Aib8 C"); axes[1,0].legend()
    axes[1,1].hist(cat_dist[half:], bins=50, alpha=0.7)
    axes[1,1].axvline(x=3.5, color='gray', linestyle='--')
    axes[1,1].set_title("Catalytic distance distribution (last half)")
    for ax in axes.flat: ax.set_xlabel("Time (ns)")
    plt.suptitle(f"{label.upper()}: {time_ns[-1]:.1f} ns trajectory", fontsize=14)
    plt.tight_layout()
    plt.savefig(f"{OUT}/{label}_quick.png", dpi=100)
    plt.close()
    print(f"  Plot: {OUT}/{label}_quick.png")

print("\nDone!")
