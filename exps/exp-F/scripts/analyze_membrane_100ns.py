#!/usr/bin/env python3
"""
Comprehensive analysis of membrane system 100 ns production MD.
Uses step=10 (579 frames) for trajectory analysis to keep runtime reasonable.

Usage:
    cd exps/exp-F && python scripts/analyze_membrane_100ns.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.rms import RMSF

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/membrane_build/system_final.prmtop"
DCD = f"{EXP_F}/md/membrane_rep1/prod_v2.dcd"
OUTDIR = f"{EXP_F}/md/membrane_rep1/analysis"
STEP = 10  # analyze every 10th frame

os.makedirs(OUTDIR, exist_ok=True)

print("Loading universe...")
u = mda.Universe(PRMTOP, DCD)
ref = u.copy()
ref.trajectory[0]
n_frames = len(u.trajectory)
n_analyzed = (n_frames + STEP - 1) // STEP
print(f"Frames: {n_frames} -> analyzing {n_analyzed} (step={STEP})")
print(f"Atoms: {len(u.atoms)}")

# Define selections
protein = u.select_atoms('protein')
peptide = u.select_atoms('resid 1114-1139')  # GLP-1 peptide
lnk = u.select_atoms('resid 1140')  # Linker
lipid = u.select_atoms('resname PA PC OL CHL')
water = u.select_atoms('resname WAT')

print(f"Protein: {len(protein)} atoms")
print(f"Peptide: {len(peptide)} atoms ({len(peptide.residues)} res)")
print(f"LNK: {len(lnk)} atoms")
print(f"Lipid: {len(lipid)} atoms")
print(f"Water: {len(water)} atoms")

# --- 1. Align trajectory on protein CA ---
ca_all = protein.select_atoms('name CA')
print(f"\nAligning on {len(ca_all)} CA atoms (step={STEP})...")
aligner = align.AlignTraj(u, ref, select='name CA and protein', in_memory=False).run(step=STEP)

# --- 2. RMSD Analysis ---
print("\n=== RMSD Analysis ===")

def compute_rmsd(label, sel, ref_sel=None):
    if ref_sel is None:
        ref_sel = sel
    r = rms.RMSD(u, ref, select=sel, groupselections=None)
    r.run(step=STEP)
    vals = r.results.rmsd[:, 2]
    n = len(vals)
    drift = np.mean(vals[-n//10:]) - np.mean(vals[:n//10]) if n >= 20 else np.mean(vals[-5:]) - np.mean(vals[:5])
    print(f"  {label:20s}: mean={np.mean(vals):.3f} Å  max={np.max(vals):.3f} Å  drift={drift:+.3f} Å")
    return vals

rmsd_total = compute_rmsd("Total Protein CA", 'name CA and protein')
rmsd_peptide = compute_rmsd("Peptide CA", 'name CA and resid 1114-1139')
rmsd_receptor = compute_rmsd("Receptor CA", 'name CA and protein and not resid 1114-1140')

# Save RMSD data
time_ns = np.arange(len(rmsd_total)) * 100.0 / len(rmsd_total)
rmsd_df = pd.DataFrame({
    'frame': np.arange(0, n_frames, STEP)[:len(rmsd_total)],
    'time_ns': time_ns,
    'total_protein': rmsd_total,
    'receptor': rmsd_receptor,
    'peptide': rmsd_peptide,
})
rmsd_df.to_csv(f"{OUTDIR}/rmsd.csv", index=False)

# Plot RMSD
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(rmsd_df['time_ns'], rmsd_df['total_protein'], 'b-', alpha=0.7, label=f'Total Protein ({np.mean(rmsd_total):.2f} Å)', linewidth=0.8)
ax.plot(rmsd_df['time_ns'], rmsd_df['receptor'], 'r-', alpha=0.7, label=f'Receptor ({np.mean(rmsd_receptor):.2f} Å)', linewidth=0.8)
ax.plot(rmsd_df['time_ns'], rmsd_df['peptide'], 'g-', alpha=0.7, label=f'Peptide ({np.mean(rmsd_peptide):.2f} Å)', linewidth=0.8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('RMSD (Å)')
ax.set_title(f'Membrane System RMSD (100 ns, step={STEP})')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsd.png", dpi=150)
print(f"  Saved: {OUTDIR}/rmsd.png")

# --- 3. RMSF Analysis ---
print("\n=== RMSF Analysis ===")

def compute_rmsf(label, sel_str):
    u2 = mda.Universe(PRMTOP, DCD)
    ref2 = u2.copy()
    ref2.trajectory[0]
    align.AlignTraj(u2, ref2, select='name CA and protein', in_memory=False).run(step=STEP)
    ca = u2.select_atoms(sel_str)
    r = RMSF(ca).run(step=STEP)
    vals = r.results.rmsf
    print(f"  {label:20s}: mean={np.mean(vals):.3f} Å  max={np.max(vals):.3f} Å")
    return vals, ca.resids, ca.resnames

rmsf_receptor, rec_resids, rec_resnames = compute_rmsf("Receptor CA", 'name CA and protein and not resid 1114-1140')
rmsf_peptide, pep_resids, pep_resnames = compute_rmsf("Peptide CA", 'name CA and resid 1114-1139')

# Plot RMSF
fig, ax = plt.subplots(figsize=(14, 5))
ax.plot(rec_resids, rmsf_receptor, 'ro-', alpha=0.5, markersize=2, linewidth=0.5, label='Receptor')
ax.plot(pep_resids, rmsf_peptide, 'go-', alpha=0.7, markersize=4, linewidth=1, label='Peptide')
ax.axvline(x=1113.5, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Residue ID')
ax.set_ylabel('RMSF (Å)')
ax.set_title(f'Membrane System RMSF per Residue (100 ns, step={STEP})')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsf.png", dpi=150)
print(f"  Saved: {OUTDIR}/rmsf.png")

# --- 4. LNK Position Tracking ---
print("\n=== LNK Position Tracking ===")
print(f"  LNK atoms: {list(lnk.atoms.names)}")

lnk_com_z = []
lnk_head_z = []
lnk_tail_z = []
membrane_com_z = []
membrane_upper_z = []
membrane_lower_z = []

# Use P atoms for membrane reference
phosphorus = u.select_atoms('name P31')
print(f"  Phosphate atoms: {len(phosphorus)}")

for ts in u.trajectory[::STEP]:
    p_z = phosphorus.positions[:, 2]
    mem_com = np.mean(p_z)
    mem_upper = np.percentile(p_z, 90)
    mem_lower = np.percentile(p_z, 10)
    
    membrane_com_z.append(mem_com)
    membrane_upper_z.append(mem_upper)
    membrane_lower_z.append(mem_lower)
    
    lnk_com_z.append(lnk.center_of_mass()[2])
    
    lnk_heavy = lnk.select_atoms('not name H*')
    lnk_pos = lnk_heavy.positions[:, 2]
    lnk_head_z.append(np.min(lnk_pos))
    lnk_tail_z.append(np.max(lnk_pos))

lnk_com_z = np.array(lnk_com_z)
lnk_head_z = np.array(lnk_head_z)
lnk_tail_z = np.array(lnk_tail_z)
membrane_com_z = np.array(membrane_com_z)
membrane_upper_z = np.array(membrane_upper_z)
membrane_lower_z = np.array(membrane_lower_z)

lnk_rel_com = lnk_com_z - membrane_com_z
lnk_tail_rel = lnk_tail_z - membrane_upper_z
lnk_head_rel = lnk_head_z - membrane_upper_z

print(f"  LNK COM relative to membrane COM: {np.mean(lnk_rel_com):.2f} ± {np.std(lnk_rel_com):.2f} Å")
print(f"  LNK tail relative to membrane upper: {np.mean(lnk_tail_rel):.2f} ± {np.std(lnk_tail_rel):.2f} Å")
print(f"  LNK head relative to membrane upper: {np.mean(lnk_head_rel):.2f} ± {np.std(lnk_head_rel):.2f} Å")

# Plot LNK position
fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

axes[0].plot(time_ns, lnk_com_z, 'purple', alpha=0.7, label='LNK COM', linewidth=0.8)
axes[0].plot(time_ns, membrane_com_z, 'gray', alpha=0.5, label='Membrane COM', linewidth=0.8)
axes[0].fill_between(time_ns, membrane_lower_z, membrane_upper_z, alpha=0.2, color='blue', label='Membrane (P 10-90%)')
axes[0].set_ylabel('z (Å)')
axes[0].set_title('LNK Position in Membrane System')
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)

axes[1].plot(time_ns, lnk_rel_com, 'purple', alpha=0.7, linewidth=0.8)
axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[1].set_ylabel('LNK COM - Mem COM (Å)')
axes[1].grid(True, alpha=0.3)

axes[2].plot(time_ns, lnk_tail_rel, 'red', alpha=0.7, label='LNK tail - mem upper', linewidth=0.8)
axes[2].plot(time_ns, lnk_head_rel, 'blue', alpha=0.7, label='LNK head - mem upper', linewidth=0.8)
axes[2].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[2].set_ylabel('Relative z (Å)')
axes[2].set_xlabel('Time (ns)')
axes[2].legend(fontsize=8)
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/lnk_position.png", dpi=150)
print(f"  Saved: {OUTDIR}/lnk_position.png")

# --- 5. Membrane Properties ---
print("\n=== Membrane Properties ===")

membrane_thickness = membrane_upper_z - membrane_lower_z
print(f"  Membrane thickness (P 90-10%): {np.mean(membrane_thickness):.2f} ± {np.std(membrane_thickness):.2f} Å")

# Area per lipid
box_xy = []
for ts in u.trajectory[::STEP]:
    box = ts.dimensions[:3]
    box_xy.append(box[0] * box[1])
box_xy = np.array(box_xy)
n_lipid_molecules = len(phosphorus.residues)
area_per_lipid = box_xy / (n_lipid_molecules / 2)
print(f"  Box area: {np.mean(box_xy):.1f} ± {np.std(box_xy):.1f} Å²")
print(f"  Area per lipid: {np.mean(area_per_lipid):.2f} ± {np.std(area_per_lipid):.2f} Å²")

fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
axes[0].plot(time_ns, membrane_thickness, 'b-', alpha=0.7, linewidth=0.8)
axes[0].axhline(y=np.mean(membrane_thickness), color='r', linestyle='--', alpha=0.5)
axes[0].set_ylabel('Thickness (Å)')
axes[0].set_title('Membrane Thickness (P 90-10%)')
axes[0].grid(True, alpha=0.3)

axes[1].plot(time_ns, area_per_lipid, 'g-', alpha=0.7, linewidth=0.8)
axes[1].axhline(y=np.mean(area_per_lipid), color='r', linestyle='--', alpha=0.5)
axes[1].set_ylabel('Area per lipid (Å²)')
axes[1].set_xlabel('Time (ns)')
axes[1].set_title('Area per Lipid')
axes[1].grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/membrane_properties.png", dpi=150)
print(f"  Saved: {OUTDIR}/membrane_properties.png")

# --- 6. Energy Analysis ---
print("\n=== Energy Analysis ===")
log_path = f"{EXP_F}/md/membrane_rep1/prod_v2.log"
if os.path.exists(log_path):
    prod = pd.read_csv(log_path, skiprows=1, header=None)
    prod.columns = ['Step', 'Time_ps', 'PE', 'KE', 'TE', 'Temp', 'Vol', 'Dens', 'Speed']
    print(f"  Frames: {len(prod)}")
    print(f"  Temp: {prod['Temp'].mean():.2f} ± {prod['Temp'].std():.2f} K")
    print(f"  Density: {prod['Dens'].mean():.4f} ± {prod['Dens'].std():.4f} g/mL")
    print(f"  Volume: {prod['Vol'].mean():.1f} ± {prod['Vol'].std():.1f} nm³")
    
    pe_first = np.mean(prod['PE'][:50])
    pe_last = np.mean(prod['PE'][-50:])
    print(f"  PE drift (first 50 vs last 50): {pe_last - pe_first:.1f} kJ/mol")
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    t = prod['Time_ps'] / 1000
    axes[0].plot(t, prod['Temp'], 'r-', alpha=0.7, linewidth=0.5)
    axes[0].axhline(y=310, color='k', linestyle='--', alpha=0.3)
    axes[0].set_ylabel('T (K)')
    axes[0].set_title('Membrane System Thermodynamics (100 ns)')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].plot(t, prod['PE']/1000, 'b-', alpha=0.7, linewidth=0.5)
    axes[1].set_ylabel('PE (MJ/mol)')
    axes[1].grid(True, alpha=0.3)
    
    axes[2].plot(t, prod['Dens'], 'purple', alpha=0.7, linewidth=0.5)
    axes[2].axhline(y=1.0, color='k', linestyle='--', alpha=0.3)
    axes[2].set_ylabel('Density (g/mL)')
    axes[2].set_xlabel('Time (ns)')
    axes[2].grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/thermodynamics.png", dpi=150)
    print(f"  Saved: {OUTDIR}/thermodynamics.png")
else:
    print(f"  Log not found: {log_path}")

# --- 7. Summary CSV ---
print("\n=== Summary ===")
summary = {
    'total_frames': n_frames,
    'analyzed_frames': len(rmsd_total),
    'total_ns': 100.0,
    'receptor_rmsd_mean_A': float(np.mean(rmsd_receptor)),
    'receptor_rmsd_max_A': float(np.max(rmsd_receptor)),
    'peptide_rmsd_mean_A': float(np.mean(rmsd_peptide)),
    'peptide_rmsd_max_A': float(np.max(rmsd_peptide)),
    'receptor_rmsf_mean_A': float(np.mean(rmsf_receptor)),
    'receptor_rmsf_max_A': float(np.max(rmsf_receptor)),
    'peptide_rmsf_mean_A': float(np.mean(rmsf_peptide)),
    'peptide_rmsf_max_A': float(np.max(rmsf_peptide)),
    'lnk_tail_membrane_dist_A': float(np.mean(lnk_tail_rel)),
    'lnk_tail_membrane_std_A': float(np.std(lnk_tail_rel)),
    'membrane_thickness_A': float(np.mean(membrane_thickness)),
    'area_per_lipid_A2': float(np.mean(area_per_lipid)),
}

summary_df = pd.DataFrame([summary])
summary_df.to_csv(f"{OUTDIR}/summary.csv", index=False)
print(f"  Saved: {OUTDIR}/summary.csv")

for k, v in summary.items():
    if isinstance(v, float):
        print(f"  {k:35s}: {v:.3f}")
    else:
        print(f"  {k:35s}: {v}")

print(f"\nAnalysis complete. Output: {OUTDIR}/")
