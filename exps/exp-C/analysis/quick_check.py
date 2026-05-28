#!/usr/bin/env python3
"""Quick 100 ns analysis: HSA+FA stability, FA3 pocket distances."""
import mdtraj as md, numpy as np, os, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt

REPO = "/home/scroll/personal/semaglutide-in-silico"
OUT = f"{REPO}/exps/exp-C/analysis"
os.makedirs(OUT, exist_ok=True)

prmtop = f"{REPO}/exps/exp-C/tleap/c18_monoacid.prmtop"

results = {}
for rep in [1, 2, 3]:
    dcd = f"{REPO}/exps/exp-C/md/c18_monoacid/rep{rep}/c18_monoacid_traj.dcd"
    if not os.path.exists(dcd): continue

    print(f"Loading rep {rep}...")
    t = md.load(dcd, top=prmtop)
    half = t.n_frames // 2
    t_eq = t[half:]  # last 50 ns for analysis
    print(f"  {t_eq.n_frames} frames (last 50 ns)")

    # Selections — avoid special chars in residue names
    hsa_ca = t.topology.select("protein and name CA")
    fah_heavy = t.topology.select("resname FAH and not element H")
    fah_o1d = t.topology.select("resname FAH and name O1D")
    fah_o2d = t.topology.select("resname FAH and name O2D")
    # FA3 key residues
    r348_nh = t.topology.select("resid 348 and (name NH1 or name NH2)")
    r485_nh = t.topology.select("resid 485 and (name NH1 or name NH2)")

    # Align on HSA
    t.superpose(t[0], atom_indices=hsa_ca)

    # HSA RMSD
    hsa_rmsd = md.rmsd(t, t[0], atom_indices=hsa_ca) * 10

    # FA heavy atom RMSD (relative to first frame)
    fah_rmsd = md.rmsd(t, t[0], atom_indices=fah_heavy) * 10

    # FA carboxyl → R348/R485 distances
    d_348 = []; d_485 = []
    if len(fah_o1d) > 0 and len(fah_o2d) > 0 and len(r348_nh) > 0 and len(r485_nh) > 0:
        for frame in range(t.n_frames):
            fah_oxy = np.concatenate([t.xyz[frame, fah_o1d], t.xyz[frame, fah_o2d]])
            r348_xyz = t.xyz[frame, r348_nh]
            r485_xyz = t.xyz[frame, r485_nh]
            d_348.append(min(np.linalg.norm(o - n) for o in fah_oxy for n in r348_xyz) * 10)
            d_485.append(min(np.linalg.norm(o - n) for o in fah_oxy for n in r485_xyz) * 10)
        d_348 = np.array(d_348); d_485 = np.array(d_485)

    results[rep] = {
        'hsa_rmsd': hsa_rmsd, 'fah_rmsd': fah_rmsd,
        'd_348': d_348, 'd_485': d_485,
        'n_frames_eq': t_eq.n_frames
    }

    print(f"  HSA RMSD (last 50ns): {hsa_rmsd[half:].mean():.1f} ± {hsa_rmsd[half:].std():.1f} Å")
    print(f"  FAH RMSD (last 50ns): {fah_rmsd[half:].mean():.1f} ± {fah_rmsd[half:].std():.1f} Å")
    if len(d_348) > 0:
        print(f"  FA→R348 (last 50ns): {d_348[half:].mean():.1f} ± {d_348[half:].std():.1f} Å")
        print(f"  FA→R485 (last 50ns): {d_485[half:].mean():.1f} ± {d_485[half:].std():.1f} Å")

# Plot
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
for i, rep in enumerate([1,2,3]):
    r = results[rep]
    time_full = np.arange(len(r['hsa_rmsd'])) * 0.1
    half = len(r['hsa_rmsd']) // 2

    axes[0,i].plot(time_full, r['hsa_rmsd'], alpha=0.7)
    axes[0,i].axvline(x=time_full[half], color='gray', linestyle='--', alpha=0.5)
    axes[0,i].set_title(f"Rep {rep}: HSA RMSD (μ={r['hsa_rmsd'][half:].mean():.1f}Å)")
    axes[0,i].set_ylabel("Å"); axes[0,i].set_xlabel("Time (ns)")

    axes[1,i].plot(time_full, r['fah_rmsd'], alpha=0.7)
    axes[1,i].axvline(x=time_full[half], color='gray', linestyle='--', alpha=0.5)
    axes[1,i].set_title(f"Rep {rep}: FA RMSD (μ={r['fah_rmsd'][half:].mean():.1f}Å)")
    axes[1,i].set_ylabel("Å"); axes[1,i].set_xlabel("Time (ns)")

plt.suptitle("exp-C: C18 monoacid + HSA — 100 ns validation (3 replicas)", fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f"{OUT}/quick_check.png", dpi=120)
plt.close()

# R348/R485 distance plot
if all(len(results[r].get('d_348',[])) > 0 for r in [1,2,3]):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for rep in [1,2,3]:
        r = results[rep]
        time_full = np.arange(len(r['d_348'])) * 0.1
        axes[0].plot(time_full, r['d_348'], alpha=0.7, label=f'Rep {rep}')
        axes[1].plot(time_full, r['d_485'], alpha=0.7, label=f'Rep {rep}')
    axes[0].set_title("FA carboxyl → R348 distance"); axes[0].set_ylabel("Å"); axes[0].legend()
    axes[1].set_title("FA carboxyl → R485 distance"); axes[1].set_ylabel("Å"); axes[1].legend()
    for ax in axes: ax.set_xlabel("Time (ns)")
    plt.suptitle("exp-C: FA3 pocket key distances", fontweight='bold')
    plt.tight_layout()
    plt.savefig(f"{OUT}/fa3_distances.png", dpi=120)
    plt.close()

print(f"\nPlots saved: {OUT}/quick_check.png, {OUT}/fa3_distances.png")
print("Done!")
