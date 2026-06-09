#!/usr/bin/env python3
"""
Analysis of OLD membrane system 108 ns production MD (with LNK defect).
Uses Merge+dummy atoms to handle topology mismatch (312501 atoms in DCD).

Metrics that are still meaningful despite LNK defect:
  - Receptor/TM dynamics (TM6 displacement, RMSD)
  - Peptide-receptor interactions
  - Membrane properties (thickness, area per lipid)
  - LNK position (rough, but note 4 extra atoms in LNK)

Metrics that are CORRUPTED by LNK defect:
  - LNK-specific internal geometry (bond angles, torsions)
  - LNK-protein interactions involving the extra atoms
  - LNK internal flexibility/RMSF
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
PDB = f"{EXP_F}/membrane_build/system_final.pdb"
DCD_OLD = f"{EXP_F}/md/membrane_rep1/old_defective/prod_v2.dcd"
OUTDIR = f"{EXP_F}/md/membrane_rep1/old_defective/analysis"
STEP = 10  # analyze every 10th frame

os.makedirs(OUTDIR, exist_ok=True)

# Build universe with dummy atoms for old DCD
print("Building universe for old DCD (312501 atoms)...")
ref = mda.Universe(PDB)
print(f"  PDB reference atoms: {ref.atoms.n_atoms}")

# 25 dummy atoms
dummy = mda.Universe.empty(n_atoms=25, n_residues=1, trajectory=False)
dummy.add_TopologyAttr('names', [f'D{i+1}' for i in range(25)])
dummy.add_TopologyAttr('types', ['D'] * 25)
dummy.add_TopologyAttr('resnames', ['DUM'])
dummy.add_TopologyAttr('resids', [99999])
dummy.add_TopologyAttr('segids', ['D'])
dummy.add_TopologyAttr('masses', [1.0] * 25)
dummy.add_TopologyAttr('charges', [0.0] * 25)

u = mda.Merge(ref.atoms, dummy.atoms)
u.load_new(DCD_OLD)
ref_u = u.copy()
ref_u.trajectory[0]

n_frames = len(u.trajectory)
n_analyzed = (n_frames + STEP - 1) // STEP
print(f"  Frames: {n_frames} -> analyzing {n_analyzed} (step={STEP})")
print(f"  Total atoms: {u.atoms.n_atoms}")

# --- Selections ---
receptor = u.select_atoms('resid 730-1113')
peptide = u.select_atoms('resid 1114-1139')
lnk = u.select_atoms('resid 1140')
phosphorus = u.select_atoms('name P31')

print(f"  Receptor: {len(receptor)} atoms ({len(receptor.residues)} residues)")
print(f"  Peptide: {len(peptide)} atoms ({len(peptide.residues)} residues)")
print(f"  LNK: {len(lnk)} atoms (NOTE: includes 4 extra garbage atoms)")
print(f"  Phosphate heads: {len(phosphorus)} atoms")

# --- 1. RMSD Analysis ---
print("\n=== RMSD Analysis ===")

def compute_rmsd(label, sel):
    r = rms.RMSD(u, ref_u, select=sel, groupselections=None)
    r.run(step=STEP)
    vals = r.results.rmsd[:, 2]
    n = len(vals)
    drift = np.mean(vals[-n//10:]) - np.mean(vals[:n//10]) if n >= 20 else np.mean(vals[-5:]) - np.mean(vals[:5])
    print(f"  {label:25s}: mean={np.mean(vals):.3f} Å  max={np.max(vals):.3f} Å  drift={drift:+.3f} Å")
    return vals

# Align on receptor CA
aligner = align.AlignTraj(u, ref_u, select='name CA and resid 730-1113', in_memory=False).run(step=STEP)

rmsd_receptor = compute_rmsd("Receptor CA", 'name CA and resid 730-1113')
rmsd_peptide = compute_rmsd("Peptide CA", 'name CA and resid 1114-1139')
rmsd_total = compute_rmsd("Total Protein CA", 'name CA and (resid 730-1139)')

# Save
actual_ns = n_frames * 10.0 / 1000  # DCD saved every 5000 steps = 10 ps
time_ns = np.arange(len(rmsd_total)) * actual_ns / len(rmsd_total)
rmsd_df = pd.DataFrame({
    'frame': np.arange(0, n_frames, STEP)[:len(rmsd_total)],
    'time_ns': time_ns,
    'total_protein': rmsd_total,
    'receptor': rmsd_receptor,
    'peptide': rmsd_peptide,
})
rmsd_df.to_csv(f"{OUTDIR}/rmsd.csv", index=False)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(rmsd_df['time_ns'], rmsd_df['total_protein'], 'b-', alpha=0.7, label=f'Total ({np.mean(rmsd_total):.2f} Å)', linewidth=0.8)
ax.plot(rmsd_df['time_ns'], rmsd_df['receptor'], 'r-', alpha=0.7, label=f'Receptor ({np.mean(rmsd_receptor):.2f} Å)', linewidth=0.8)
ax.plot(rmsd_df['time_ns'], rmsd_df['peptide'], 'g-', alpha=0.7, label=f'Peptide ({np.mean(rmsd_peptide):.2f} Å)', linewidth=0.8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('RMSD (Å)')
ax.set_title(f'OLD Membrane System RMSD (108 ns, LNK-defective)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsd.png", dpi=150)
print(f"  Saved: {OUTDIR}/rmsd.png")

# --- 2. RMSF ---
print("\n=== RMSF Analysis ===")

def compute_rmsf(label, sel_str):
    u2 = mda.Merge(ref.atoms, dummy.atoms)
    u2.load_new(DCD_OLD)
    ref2 = u2.copy()
    ref2.trajectory[0]
    align.AlignTraj(u2, ref2, select='name CA and resid 730-1113', in_memory=False).run(step=STEP)
    ca = u2.select_atoms(sel_str)
    r = RMSF(ca).run(step=STEP)
    vals = r.results.rmsf
    print(f"  {label:25s}: mean={np.mean(vals):.3f} Å  max={np.max(vals):.3f} Å")
    return vals, ca.resids

rmsf_receptor, rec_resids = compute_rmsf("Receptor CA", 'name CA and resid 730-1113')
rmsf_peptide, pep_resids = compute_rmsf("Peptide CA", 'name CA and resid 1114-1139')

fig, ax = plt.subplots(figsize=(14, 5))
ax.plot(rec_resids, rmsf_receptor, 'ro-', alpha=0.5, markersize=2, linewidth=0.5, label='Receptor')
ax.plot(pep_resids, rmsf_peptide, 'go-', alpha=0.7, markersize=4, linewidth=1, label='Peptide')
ax.axvline(x=1113.5, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Residue ID')
ax.set_ylabel('RMSF (Å)')
ax.set_title(f'OLD Membrane System RMSF (108 ns, LNK-defective)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsf.png", dpi=150)
print(f"  Saved: {OUTDIR}/rmsf.png")

# --- 3. TM6 Displacement ---
print("\n=== TM6 Displacement ===")
# TM6 = resid 1079-1113
tm6 = u.select_atoms('resid 1079-1113 and name CA')
tm1 = u.select_atoms('resid 841-886 and name CA')
tm2 = u.select_atoms('resid 933-972 and name CA')
tm3 = u.select_atoms('resid 1014-1051 and name CA')

print(f"  TM6 CA atoms: {len(tm6)}")
print(f"  TM1 CA atoms: {len(tm1)}")
print(f"  TM2 CA atoms: {len(tm2)}")
print(f"  TM3 CA atoms: {len(tm3)}")

tm6_z = []
tm1_z = []
tm2_z = []
tm3_z = []
membrane_z = []

for ts in u.trajectory[::STEP]:
    tm6_z.append(np.mean(tm6.positions[:, 2]))
    tm1_z.append(np.mean(tm1.positions[:, 2]))
    tm2_z.append(np.mean(tm2.positions[:, 2]))
    tm3_z.append(np.mean(tm3.positions[:, 2]))
    membrane_z.append(np.mean(phosphorus.positions[:, 2]))

tm6_z = np.array(tm6_z)
tm1_z = np.array(tm1_z)
tm2_z = np.array(tm2_z)
tm3_z = np.array(tm3_z)
membrane_z = np.array(membrane_z)

tm6_rel = tm6_z - membrane_z
tm1_rel = tm1_z - membrane_z
tm2_rel = tm2_z - membrane_z
tm3_rel = tm3_z - membrane_z

print(f"  TM6 z relative to membrane: {np.mean(tm6_rel):.2f} ± {np.std(tm6_rel):.2f} Å")
print(f"  TM1 z relative to membrane: {np.mean(tm1_rel):.2f} ± {np.std(tm1_rel):.2f} Å")
print(f"  TM2 z relative to membrane: {np.mean(tm2_rel):.2f} ± {np.std(tm2_rel):.2f} Å")
print(f"  TM3 z relative to membrane: {np.mean(tm3_rel):.2f} ± {np.std(tm3_rel):.2f} Å")

# TM6 displacement = distance from average position of TM1+TM2+TM3
tm_ref = (tm1_z + tm2_z + tm3_z) / 3
tm6_displacement = tm6_z - tm_ref
print(f"  TM6 displacement (from TM1/2/3 avg): {np.mean(tm6_displacement):.2f} ± {np.std(tm6_displacement):.2f} Å")

fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
axes[0].plot(time_ns, tm6_rel, 'purple', alpha=0.7, label='TM6', linewidth=0.8)
axes[0].plot(time_ns, tm1_rel, 'blue', alpha=0.5, label='TM1', linewidth=0.8)
axes[0].plot(time_ns, tm2_rel, 'green', alpha=0.5, label='TM2', linewidth=0.8)
axes[0].plot(time_ns, tm3_rel, 'orange', alpha=0.5, label='TM3', linewidth=0.8)
axes[0].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[0].set_ylabel('z - Membrane COM (Å)')
axes[0].set_title('TM Helix Positions Relative to Membrane')
axes[0].legend(fontsize=8)
axes[0].grid(True, alpha=0.3)

axes[1].plot(time_ns, tm6_displacement, 'purple', alpha=0.7, linewidth=0.8)
axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[1].set_ylabel('TM6 - (TM1+TM2+TM3)/3 (Å)')
axes[1].set_xlabel('Time (ns)')
axes[1].set_title('TM6 Displacement (Outward Movement Indicator)')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/tm6_displacement.png", dpi=150)
print(f"  Saved: {OUTDIR}/tm6_displacement.png")

# --- 4. LNK Position (rough, note defect) ---
print("\n=== LNK Position (rough, LNK-defective) ===")
lnk_com_z = []
lnk_tail_z = []
lnk_head_z = []
membrane_upper_z = []
membrane_lower_z = []

for ts in u.trajectory[::STEP]:
    p_z = phosphorus.positions[:, 2]
    mem_upper = np.percentile(p_z, 90)
    mem_lower = np.percentile(p_z, 10)
    mem_com = np.mean(p_z)
    
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
membrane_upper_z = np.array(membrane_upper_z)
membrane_lower_z = np.array(membrane_lower_z)

lnk_rel_com = lnk_com_z - membrane_z
lnk_tail_rel = lnk_tail_z - membrane_upper_z
lnk_head_rel = lnk_head_z - membrane_upper_z

print(f"  LNK COM relative to membrane COM: {np.mean(lnk_rel_com):.2f} ± {np.std(lnk_rel_com):.2f} Å")
print(f"  LNK tail relative to membrane upper: {np.mean(lnk_tail_rel):.2f} ± {np.std(lnk_tail_rel):.2f} Å")
print(f"  LNK head relative to membrane upper: {np.mean(lnk_head_rel):.2f} ± {np.std(lnk_head_rel):.2f} Å")

fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
axes[0].plot(time_ns, lnk_com_z, 'purple', alpha=0.7, label='LNK COM', linewidth=0.8)
axes[0].plot(time_ns, membrane_z, 'gray', alpha=0.5, label='Membrane COM', linewidth=0.8)
axes[0].fill_between(time_ns, membrane_lower_z, membrane_upper_z, alpha=0.2, color='blue', label='Membrane (P 10-90%)')
axes[0].set_ylabel('z (Å)')
axes[0].set_title('LNK Position (108 ns, LNK-defective: 4 extra atoms)')
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

# --- 5. Fatty Acid Chain Angle ---
print("\n=== Fatty Acid Chain Angle ===")
# LNK atoms: C47-C48-C49-C50-C51-C52-C53-C54-C55 (last 9 carbons of chain)
# Measure angle between vector from C47->C55 and membrane normal (z-axis)
# Note: in defective LNK, there are extra atoms, so we use named atoms

try:
    c47 = u.select_atoms('resid 1140 and name C47')
    c55 = u.select_atoms('resid 1140 and name C55')
    
    if len(c47) == 1 and len(c55) == 1:
        angles = []
        for ts in u.trajectory[::STEP]:
            v = c55.positions[0] - c47.positions[0]
            v_norm = v / np.linalg.norm(v)
            # Angle with z-axis
            cos_theta = abs(v_norm[2])  # abs because chain can point either way
            angle_deg = np.degrees(np.arccos(np.clip(cos_theta, -1, 1)))
            angles.append(angle_deg)
        
        angles = np.array(angles)
        print(f"  Fatty acid chain angle (C47-C55 vs z-axis): {np.mean(angles):.1f} ± {np.std(angles):.1f}°")
        print(f"  Range: {np.min(angles):.1f}° - {np.max(angles):.1f}°")
        
        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(time_ns, angles, 'purple', alpha=0.7, linewidth=0.8)
        ax.axhline(y=np.mean(angles), color='r', linestyle='--', alpha=0.5)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Angle (°)')
        ax.set_title('Fatty Acid Chain Angle (C47-C55 vs Membrane Normal)')
        ax.set_ylim(0, 90)
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{OUTDIR}/fatty_acid_angle.png", dpi=150)
        print(f"  Saved: {OUTDIR}/fatty_acid_angle.png")
    else:
        print(f"  Could not find C47/C55 (found {len(c47)}/{len(c55)})")
except Exception as e:
    print(f"  Error computing fatty acid angle: {e}")

# --- 6. Membrane Properties ---
print("\n=== Membrane Properties ===")
membrane_thickness = membrane_upper_z - membrane_lower_z
print(f"  Membrane thickness (P 90-10%): {np.mean(membrane_thickness):.2f} ± {np.std(membrane_thickness):.2f} Å")

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

# --- 7. Energy ---
print("\n=== Energy Analysis ===")
log_path = f"{EXP_F}/md/membrane_rep1/old_defective/prod_v2.log"
if os.path.exists(log_path):
    prod = pd.read_csv(log_path, sep=',', skiprows=1, header=None, engine='python')
    prod.columns = ['Step', 'Time_ps', 'PE', 'KE', 'TE', 'Temp', 'Vol', 'Dens', 'Speed']
    # Convert to numeric and drop any malformed rows
    for col in prod.columns:
        prod[col] = pd.to_numeric(prod[col], errors='coerce')
    prod = prod.dropna()
    print(f"  Frames: {len(prod)}")
    print(f"  Temp: {prod['Temp'].mean():.2f} ± {prod['Temp'].std():.2f} K")
    print(f"  Density: {prod['Dens'].mean():.4f} ± {prod['Dens'].std():.4f} g/mL")
    
    pe_first = np.mean(prod['PE'][:50])
    pe_last = np.mean(prod['PE'][-50:])
    pe_drift = pe_last - pe_first
    print(f"  PE drift (first 50 vs last 50): {pe_drift:.1f} kJ/mol")
    
    # PE jump analysis
    pe = prod['PE'].values
    jumps = np.diff(pe)
    big_jumps = np.sum(np.abs(jumps) > 1000)
    print(f"  PE jumps >1000 kJ/mol: {big_jumps} ({100*big_jumps/len(jumps):.1f}%)")
    
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    t = prod['Time_ps'] / 1000
    axes[0].plot(t, prod['Temp'], 'r-', alpha=0.7, linewidth=0.5)
    axes[0].axhline(y=310, color='k', linestyle='--', alpha=0.3)
    axes[0].set_ylabel('T (K)')
    axes[0].set_title('OLD Membrane Thermodynamics (108 ns, LNK-defective)')
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

# --- 8. Summary ---
print("\n=== Summary ===")
summary = {
    'note': 'OLD_DATA_LNK_DEFECTIVE',
    'total_frames': n_frames,
    'actual_ns': round(actual_ns, 2),
    'receptor_rmsd_mean_A': float(np.mean(rmsd_receptor)),
    'receptor_rmsd_max_A': float(np.max(rmsd_receptor)),
    'peptide_rmsd_mean_A': float(np.mean(rmsd_peptide)),
    'peptide_rmsd_max_A': float(np.max(rmsd_peptide)),
    'receptor_rmsf_mean_A': float(np.mean(rmsf_receptor)),
    'receptor_rmsf_max_A': float(np.max(rmsf_receptor)),
    'peptide_rmsf_mean_A': float(np.mean(rmsf_peptide)),
    'peptide_rmsf_max_A': float(np.max(rmsf_peptide)),
    'tm6_displacement_mean_A': float(np.mean(tm6_displacement)),
    'tm6_displacement_std_A': float(np.std(tm6_displacement)),
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
print("WARNING: This data has LNK defects. LNK-specific metrics are unreliable.")
