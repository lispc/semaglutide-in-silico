#!/usr/bin/env python3
"""
Corrected membrane analysis for LNK tail z-position.
Memory-efficient: loads only needed atoms.
"""
import sys, os
import numpy as np
import mdtraj as md

REPO = "/home/scroll/personal/semaglutide-in-silico"
PRMTOP = f"{REPO}/exps/exp-F/membrane_build/system_final.prmtop"
DCD = f"{REPO}/exps/exp-F/md/membrane_rep1/prod_v2.dcd"
OUT = f"{REPO}/exps/exp-F/analysis"
os.makedirs(OUT, exist_ok=True)

print("Loading topology...")
top = md.load_prmtop(PRMTOP)

# === Atom selections ===
# Lipid atoms (residues >= 1140)
lipid_atoms = [a.index for a in top.atoms if a.residue.index >= 1140]
print(f"  Lipid atoms: {len(lipid_atoms)}")

# LNK residue
lnk_atoms = [a.index for a in top.atoms if a.residue.index == 1139]
print(f"  LNK atoms: {len(lnk_atoms)}")

# LNK tail atoms (diacid end)
tail_names = ['C55', 'O56', 'O57', 'C58', 'O59']
tail_atoms = [a.index for a in top.atoms if a.residue.index == 1139 and a.name in tail_names]
print(f"  LNK tail atoms: {len(tail_atoms)}")
for a in top.atoms:
    if a.residue.index == 1139 and a.name in tail_names:
        print(f"    {a.name}: mass={a.element.mass:.2f}")

# Peptide CA
peptide_ca = [a.index for a in top.atoms if a.name == 'CA' and 1113 <= a.residue.index <= 1138]
print(f"  Peptide CA: {len(peptide_ca)}")

# Receptor CA
receptor_ca = [a.index for a in top.atoms if a.name == 'CA' and a.residue.index < 1113]
print(f"  Receptor CA: {len(receptor_ca)}")

# Lipid P atoms for leaflet identification
p_atoms = [a.index for a in top.atoms if a.element.symbol == 'P' and a.residue.index >= 1140]
print(f"  Lipid P atoms: {len(p_atoms)}")

# Combine all needed atoms for loading
all_needed = list(set(lipid_atoms + lnk_atoms + peptide_ca + receptor_ca + p_atoms))
print(f"  Total atoms to load: {len(all_needed)}")

# Build index mapping from original to loaded
loaded_map = {orig: i for i, orig in enumerate(all_needed)}

# === Load trajectory (subset of atoms) ===
print("\nLoading trajectory...")
traj = md.load_dcd(DCD, top=PRMTOP, atom_indices=all_needed)
print(f"  Frames: {traj.n_frames}")
print(f"  Box (first): {traj.unitcell_lengths[0]}")

# Make molecules whole
# traj = traj.image_molecules(inplace=False)
print("  Using raw DCD coordinates (OpenMM default: unwrapped)")

# Map indices
lipid_idx = [loaded_map[i] for i in lipid_atoms]
lnk_idx = [loaded_map[i] for i in lnk_atoms]
tail_idx = [loaded_map[i] for i in tail_atoms]
pep_idx = [loaded_map[i] for i in peptide_ca]
rec_idx = [loaded_map[i] for i in receptor_ca]
p_idx = [loaded_map[i] for i in p_atoms]

# Pre-compute masses
lipid_masses = np.array([top.atom(i).element.mass for i in lipid_atoms])
lnk_masses = np.array([top.atom(i).element.mass for i in lnk_atoms])
tail_masses = np.array([top.atom(i).element.mass for i in tail_atoms])
pep_masses = np.array([top.atom(i).element.mass for i in peptide_ca])
rec_masses = np.array([top.atom(i).element.mass for i in receptor_ca])

n_frames = traj.n_frames
membrane_com_z = np.zeros(n_frames)
lnk_com_z = np.zeros(n_frames)
lnk_tail_com_z = np.zeros(n_frames)
peptide_com_z = np.zeros(n_frames)
receptor_com_z = np.zeros(n_frames)

print("\nComputing per-frame COM...")
for i in range(n_frames):
    coords = traj.xyz[i]
    membrane_com_z[i] = np.average(coords[lipid_idx, 2], weights=lipid_masses)
    lnk_com_z[i] = np.average(coords[lnk_idx, 2], weights=lnk_masses)
    lnk_tail_com_z[i] = np.average(coords[tail_idx, 2], weights=tail_masses)
    peptide_com_z[i] = np.average(coords[pep_idx, 2], weights=pep_masses)
    receptor_com_z[i] = np.average(coords[rec_idx, 2], weights=rec_masses)
    if (i + 1) % 100 == 0:
        print(f"  Frame {i+1}/{n_frames}", end='\r', flush=True)
print(f"\n  Done")

# === Relative positions ===
lnk_rel = lnk_com_z - membrane_com_z
lnk_tail_rel = lnk_tail_com_z - membrane_com_z
peptide_rel = peptide_com_z - membrane_com_z
receptor_rel = receptor_com_z - membrane_com_z

# === Bilayer thickness and surfaces ===
p_z_first = traj.xyz[0, p_idx, 2]
mid_z = np.median(p_z_first)
upper_mask = p_z_first > mid_z
lower_mask = p_z_first <= mid_z

thickness = np.zeros(n_frames)
upper_surface = np.zeros(n_frames)
lower_surface = np.zeros(n_frames)

for i in range(n_frames):
    z_rel = traj.xyz[i, p_idx, 2] - membrane_com_z[i]
    thickness[i] = np.mean(z_rel[upper_mask]) - np.mean(z_rel[lower_mask])
    upper_surface[i] = np.mean(z_rel[upper_mask])
    lower_surface[i] = np.mean(z_rel[lower_mask])

# === Summary ===
print("\n" + "="*60)
print("CORRECTED MEMBRANE ANALYSIS")
print("="*60)

print(f"\n[Absolute positions (raw DCD coordinates)]")
print(f"  Membrane COM z:    {np.mean(membrane_com_z)*10:.2f} ± {np.std(membrane_com_z)*10:.2f} Å")
print(f"  Receptor CA COM z: {np.mean(receptor_com_z)*10:.2f} ± {np.std(receptor_com_z)*10:.2f} Å")
print(f"  Peptide CA COM z:  {np.mean(peptide_com_z)*10:.2f} ± {np.std(peptide_com_z)*10:.2f} Å")
print(f"  LNK COM z:         {np.mean(lnk_com_z)*10:.2f} ± {np.std(lnk_com_z)*10:.2f} Å")
print(f"  LNK tail COM z:    {np.mean(lnk_tail_com_z)*10:.2f} ± {np.std(lnk_tail_com_z)*10:.2f} Å")

print(f"\n[Relative to membrane COM (corrected)]")
print(f"  Receptor CA:       {np.mean(receptor_rel)*10:.2f} ± {np.std(receptor_rel)*10:.2f} Å")
print(f"  Peptide CA:        {np.mean(peptide_rel)*10:.2f} ± {np.std(peptide_rel)*10:.2f} Å")
print(f"  LNK COM:           {np.mean(lnk_rel)*10:.2f} ± {np.std(lnk_rel)*10:.2f} Å")
print(f"  LNK tail (diacid): {np.mean(lnk_tail_rel)*10:.2f} ± {np.std(lnk_tail_rel)*10:.2f} Å")

print(f"\n[Membrane geometry]")
print(f"  Bilayer thickness (P-P): {np.mean(thickness)*10:.2f} ± {np.std(thickness)*10:.2f} Å")
print(f"  Upper leaflet surface:   {np.mean(upper_surface)*10:.2f} ± {np.std(upper_surface)*10:.2f} Å")
print(f"  Lower leaflet surface:   {np.mean(lower_surface)*10:.2f} ± {np.std(lower_surface)*10:.2f} Å")

print(f"\n[LNK tail distance from membrane]")
dist_to_upper = lnk_tail_rel - upper_surface
print(f"  Distance to upper surface: {np.mean(dist_to_upper)*10:.2f} ± {np.std(dist_to_upper)*10:.2f} Å")

if np.mean(dist_to_upper) > 0:
    print(f"  → LNK tail is {np.mean(dist_to_upper)*10:.1f} Å ABOVE the upper membrane surface")
else:
    print(f"  → LNK tail is {np.mean(np.abs(dist_to_upper))*10:.1f} Å BELOW/BURIED in membrane")

# === Save data ===
time_ns = np.arange(n_frames) * 0.01  # 10 ps/frame

np.savetxt(f"{OUT}/membrane_corrected_com_z.dat",
           np.column_stack([time_ns, membrane_com_z*10, receptor_com_z*10, peptide_com_z*10, lnk_com_z*10, lnk_tail_com_z*10]),
           header='Time(ns) MembraneCOM(A) ReceptorCOM(A) PeptideCOM(A) LNKCOM(A) LNKTailCOM(A)')

np.savetxt(f"{OUT}/membrane_corrected_rel_z.dat",
           np.column_stack([time_ns, receptor_rel*10, peptide_rel*10, lnk_rel*10, lnk_tail_rel*10, upper_surface*10, lower_surface*10]),
           header='Time(ns) ReceptorRel(A) PeptideRel(A) LNKRel(A) LNKTailRel(A) UpperSurface(A) LowerSurface(A)')

# === Plot ===
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

axes[0].plot(time_ns, membrane_com_z*10, label='Membrane COM', color='gray', alpha=0.7)
axes[0].plot(time_ns, receptor_com_z*10, label='Receptor CA COM', color='tab:blue', alpha=0.7)
axes[0].plot(time_ns, peptide_com_z*10, label='Peptide CA COM', color='tab:red', alpha=0.7)
axes[0].plot(time_ns, lnk_tail_com_z*10, label='LNK tail COM', color='tab:orange', alpha=0.7)
axes[0].set_ylabel('z-position (Å)')
axes[0].set_title('Absolute z-positions (raw DCD) - showing ~10 nm drift!')
axes[0].legend(loc='upper left', fontsize=8)
axes[0].grid(True, alpha=0.3)

axes[1].plot(time_ns, receptor_rel*10, label='Receptor CA COM', color='tab:blue', alpha=0.7)
axes[1].plot(time_ns, peptide_rel*10, label='Peptide CA COM', color='tab:red', alpha=0.7)
axes[1].plot(time_ns, lnk_rel*10, label='LNK COM', color='tab:green', alpha=0.7)
axes[1].plot(time_ns, lnk_tail_rel*10, label='LNK tail COM', color='tab:orange', alpha=0.7)
axes[1].axhline(y=0, color='gray', linestyle='--', alpha=0.5)
axes[1].fill_between(time_ns, np.mean(lower_surface)*10, np.mean(upper_surface)*10,
                      color='gray', alpha=0.2, label='Membrane')
axes[1].set_ylabel('z relative to membrane COM (Å)')
axes[1].set_title('Corrected positions (membrane COM recentered)')
axes[1].legend(loc='upper left', fontsize=8)
axes[1].grid(True, alpha=0.3)

axes[2].plot(time_ns, lnk_tail_rel*10, color='tab:orange', linewidth=0.8)
axes[2].axhline(y=np.mean(upper_surface)*10, color='blue', linestyle='--', alpha=0.5,
                label=f'Upper surface ({np.mean(upper_surface)*10:.1f} Å)')
axes[2].axhline(y=0, color='gray', linestyle='--', alpha=0.3, label='Membrane center')
axes[2].axhline(y=np.mean(lower_surface)*10, color='blue', linestyle=':', alpha=0.5,
                label=f'Lower surface ({np.mean(lower_surface)*10:.1f} Å)')
axes[2].fill_between(time_ns, np.mean(lower_surface)*10, np.mean(upper_surface)*10, color='gray', alpha=0.15)
axes[2].set_ylabel('LNK tail z rel to membrane COM (Å)')
axes[2].set_xlabel('Time (ns)')
axes[2].set_title(f'LNK tail position (corrected): {np.mean(dist_to_upper)*10:.1f} Å above membrane')
axes[2].legend(loc='upper right', fontsize=8)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUT}/membrane_corrected_lnk_z.png", dpi=150)
plt.close()

print(f"\nPlot saved: {OUT}/membrane_corrected_lnk_z.png")
print("Done!")
