#!/usr/bin/env python3
"""
Full 100 ns analysis of ECD v2 production MD.

Corrected residue numbering based on system_ecd_v2.prmtop:
- ECD (GLP-1R): resid 1-100
- Peptide (semaglutide): resid 101-128
- Linker (LNK): resid 129
- HSA: resid 130-711

Lys20-LNK bond exists (k=427, req=1.38A). Key insight: ECD/Peptide/LNK form a rigid
complex connected by flexible LNK to HSA. Internal RMSD (superposition=True) measures
conformational stability; COM distances measure domain motion via linker flexibility.

Usage:
    cd exps/exp-F/md/ecd_v2 && python ../../scripts/analyze_ecd_v2_100ns.py
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
PRMTOP = f"{EXP_F}/minimal_model/system_ecd_v2.prmtop"
DCD = f"{EXP_F}/md/ecd_v2/prod.dcd"
OUTDIR = f"{EXP_F}/md/ecd_v2/analysis_100ns"
STEP = 5

os.makedirs(OUTDIR, exist_ok=True)

print("Loading universe...")
u = mda.Universe(PRMTOP, DCD)
ref = u.copy()
ref.trajectory[0]
n_frames = len(u.trajectory)
n_analyzed = len(u.trajectory[::STEP])
print(f"Frames: {n_frames} -> analyzing {n_analyzed} (step={STEP})")

# Components
ecd = u.select_atoms('resid 1-100')
peptide = u.select_atoms('resid 101-128')
linker = u.select_atoms('resid 129')
hsa = u.select_atoms('resid 130-711')

print(f"ECD: {len(ecd)} atoms")
print(f"Peptide: {len(peptide)} atoms")
print(f"Linker: {len(linker)} atoms")
print(f"HSA: {len(hsa)} atoms")

# === 0. Verify Lys20-LNK bond ===
print("\n=== Bond Verification ===")
lys_nz = u.select_atoms('resid 117 and name NZ')
lnk_c = u.select_atoms('resid 129 and name C')
u.trajectory[0]
d0 = np.linalg.norm(lys_nz.positions[0] - lnk_c.positions[0])
u.trajectory[-1]
d1 = np.linalg.norm(lys_nz.positions[0] - lnk_c.positions[0])
print(f"  Lys20 NZ - LNK C distance: frame0={d0:.3f}Å, last={d1:.3f}Å")
if d0 < 2.0 and d1 < 2.0:
    print("  ✅ Bond intact (amide bond ~1.4Å)")
else:
    print("  ⚠️ Bond may be broken!")

# --- 1. Internal RMSD (superposition=True, measures conformational stability) ---
print("\n=== Internal RMSD (conformational stability) ===")

def internal_rmsd(label, sel_str):
    u2 = mda.Universe(PRMTOP, DCD)
    ref2 = u2.copy()
    ref2.trajectory[0]
    sel = u2.select_atoms(sel_str)
    ref_sel = ref2.select_atoms(sel_str)
    
    vals = []
    for ts in u2.trajectory[::STEP]:
        r = rms.rmsd(sel.positions, ref_sel.positions, superposition=True)
        vals.append(r)
    vals = np.array(vals)
    n = len(vals)
    drift = np.mean(vals[-n//10:]) - np.mean(vals[:n//10])
    print(f"  {label:12s}: mean={np.mean(vals):.3f} Å  max={np.max(vals):.3f} Å  drift={drift:+.3f} Å")
    return vals

rmsd_ecd_int = internal_rmsd("ECD", 'name CA and resid 1-100')
rmsd_pep_int = internal_rmsd("Peptide", 'name CA and resid 101-128')
rmsd_lnk_int = internal_rmsd("Linker", 'resid 129 and not name H*')
rmsd_hsa_int = internal_rmsd("HSA", 'name CA and resid 130-711')

time_ns = np.arange(len(rmsd_ecd_int)) * (n_frames * 10 / 1000) / len(rmsd_ecd_int)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(time_ns, rmsd_ecd_int, 'b-', alpha=0.7, label=f'ECD ({np.mean(rmsd_ecd_int):.2f} Å)', linewidth=0.8)
ax.plot(time_ns, rmsd_pep_int, 'g-', alpha=0.7, label=f'Peptide ({np.mean(rmsd_pep_int):.2f} Å)', linewidth=0.8)
ax.plot(time_ns, rmsd_lnk_int, 'purple', alpha=0.7, label=f'Linker ({np.mean(rmsd_lnk_int):.2f} Å)', linewidth=0.8)
ax.plot(time_ns, rmsd_hsa_int, 'r-', alpha=0.7, label=f'HSA ({np.mean(rmsd_hsa_int):.2f} Å)', linewidth=0.8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Internal RMSD (Å)')
ax.set_title('ECD v2 Internal Conformational Stability (100 ns)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsd_internal_corrected.png", dpi=150)
print(f"  Saved: {OUTDIR}/rmsd_internal_corrected.png")

# --- 2. Domain COM Distances (measures linker flexibility) ---
print("\n=== Domain COM Distances (linker flexibility) ===")

ecd_ca = u.select_atoms('name CA and resid 1-100')
pep_ca = u.select_atoms('name CA and resid 101-128')
lnk_all = u.select_atoms('resid 129')
hsa_ca = u.select_atoms('name CA and resid 130-711')

ecd_hsa_dist = []
pep_hsa_dist = []
lnk_hsa_dist = []
ecd_pep_dist = []
pep_lnk_dist = []

for ts in u.trajectory[::STEP]:
    ecd_hsa_dist.append(np.linalg.norm(ecd_ca.center_of_mass() - hsa_ca.center_of_mass()))
    pep_hsa_dist.append(np.linalg.norm(pep_ca.center_of_mass() - hsa_ca.center_of_mass()))
    lnk_hsa_dist.append(np.linalg.norm(lnk_all.center_of_mass() - hsa_ca.center_of_mass()))
    ecd_pep_dist.append(np.linalg.norm(ecd_ca.center_of_mass() - pep_ca.center_of_mass()))
    pep_lnk_dist.append(np.linalg.norm(pep_ca.center_of_mass() - lnk_all.center_of_mass()))

ecd_hsa_dist = np.array(ecd_hsa_dist)
pep_hsa_dist = np.array(pep_hsa_dist)
lnk_hsa_dist = np.array(lnk_hsa_dist)
ecd_pep_dist = np.array(ecd_pep_dist)
pep_lnk_dist = np.array(pep_lnk_dist)

print(f"  ECD-HSA:   {np.mean(ecd_hsa_dist):.1f} ± {np.std(ecd_hsa_dist):.1f} Å")
print(f"  Pep-HSA:   {np.mean(pep_hsa_dist):.1f} ± {np.std(pep_hsa_dist):.1f} Å")
print(f"  LNK-HSA:   {np.mean(lnk_hsa_dist):.1f} ± {np.std(lnk_hsa_dist):.1f} Å")
print(f"  ECD-Pep:   {np.mean(ecd_pep_dist):.1f} ± {np.std(ecd_pep_dist):.1f} Å")
print(f"  Pep-LNK:   {np.mean(pep_lnk_dist):.1f} ± {np.std(pep_lnk_dist):.1f} Å")

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(time_ns, ecd_hsa_dist, 'b-', alpha=0.7, label=f'ECD-HSA ({np.mean(ecd_hsa_dist):.1f} Å)', linewidth=0.8)
ax.plot(time_ns, pep_hsa_dist, 'g-', alpha=0.7, label=f'Pep-HSA ({np.mean(pep_hsa_dist):.1f} Å)', linewidth=0.8)
ax.plot(time_ns, lnk_hsa_dist, 'purple', alpha=0.7, label=f'LNK-HSA ({np.mean(lnk_hsa_dist):.1f} Å)', linewidth=0.8)
ax.plot(time_ns, ecd_pep_dist, 'c-', alpha=0.7, label=f'ECD-Pep ({np.mean(ecd_pep_dist):.1f} Å)', linewidth=0.8)
ax.plot(time_ns, pep_lnk_dist, 'orange', alpha=0.7, label=f'Pep-LNK ({np.mean(pep_lnk_dist):.1f} Å)', linewidth=0.8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('COM Distance (Å)')
ax.set_title('ECD v2 Domain COM Distances (100 ns)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/com_distances.png", dpi=150)
print(f"  Saved: {OUTDIR}/com_distances.png")

# Save distance data
pd.DataFrame({
    'time_ns': time_ns,
    'ecd_hsa': ecd_hsa_dist,
    'pep_hsa': pep_hsa_dist,
    'lnk_hsa': lnk_hsa_dist,
    'ecd_pep': ecd_pep_dist,
    'pep_lnk': pep_lnk_dist,
}).to_csv(f"{OUTDIR}/com_distances.csv", index=False)

# --- 3. RMSD Plateau Check ---
print("\n=== RMSD Plateau Analysis ===")

def check_plateau(vals, window=50):
    n = len(vals)
    if n < 2*window:
        window = n // 4
    early = np.mean(vals[:window])
    mid = np.mean(vals[n//2-window//2:n//2+window//2])
    late = np.mean(vals[-window:])
    return early, mid, late

for label, vals in [("ECD", rmsd_ecd_int), ("Peptide", rmsd_pep_int), ("Linker", rmsd_lnk_int), ("HSA", rmsd_hsa_int)]:
    e, m, l = check_plateau(vals)
    drift_total = l - e
    drift_mid_late = l - m
    print(f"  {label:12s}: early={e:.2f}  mid={m:.2f}  late={l:.2f}  drift={drift_total:+.3f} Å")
    if abs(drift_total) < 1.0:
        print(f"               -> PLATEAU (drift < 1.0 Å)")
    elif abs(drift_total) < 2.0:
        print(f"               -> NEAR PLATEAU (drift < 2.0 Å)")
    else:
        print(f"               -> STILL DRIFTING")

# --- 4. Internal RMSF ---
print("\n=== Internal RMSF ===")

# CRITICAL: ECD and Peptide are NOT covalently linked in this model.
# They are held by non-covalent interactions and can drift apart.
# Therefore RMSF MUST be computed with component-wise alignment:
# - ECD aligned on ECD CA only
# - Peptide aligned on Peptide CA only  
# - HSA aligned on HSA CA only
# Using any global or combined alignment artifactually inflates RMSF.

def compute_rmsf_component(label, resid_range, align_sel, rmsf_sel, heavy_only=False):
    """Load trajectory, align on specified atoms, compute RMSF."""
    print(f"  Aligning {label}...")
    u = mda.Universe(PRMTOP, DCD)
    ref = u.copy()
    ref.trajectory[0]
    align.AlignTraj(u, ref, select=align_sel, in_memory=True).run(step=STEP)
    atoms = u.select_atoms(rmsf_sel)
    r = RMSF(atoms).run(step=STEP)
    vals = r.results.rmsf
    print(f"  {label:12s}: mean={np.mean(vals):.3f} Å  max={np.max(vals):.3f} Å")
    return vals, atoms.resids, atoms.resnames

# ECD RMSF (align on ECD CA only)
rmsf_ecd, ecd_resids, ecd_resnames = compute_rmsf_component(
    "ECD", "1-100", "name CA and resid 1-100", "name CA and resid 1-100"
)

# Peptide RMSF (align on Peptide CA only)
rmsf_pep, pep_resids, pep_resnames = compute_rmsf_component(
    "Peptide", "101-128", "name CA and resid 101-128", "name CA and resid 101-128"
)

# HSA RMSF (align on HSA CA only)
rmsf_hsa, hsa_resids, hsa_resnames = compute_rmsf_component(
    "HSA", "130-711", "name CA and resid 130-711", "name CA and resid 130-711"
)

# LNK RMSF (align on LNK heavy atoms)
rmsf_lnk, lnk_resids, lnk_resnames = compute_rmsf_component(
    "Linker", "129", "resid 129 and not name H*", "resid 129 and not name H*"
)
print("\n=== HSA High RMSF Residues (>5 Å) ===")
for i, (rid, rnm, val) in enumerate(zip(hsa_resids, hsa_resnames, rmsf_hsa)):
    if val > 5.0:
        print(f"  Residue {rid} {rnm}: {val:.2f} Å")

# Top 10
print("\n=== HSA Top 10 Highest RMSF ===")
top10 = np.argsort(rmsf_hsa)[-10:][::-1]
for idx in top10:
    print(f"  Residue {hsa_resids[idx]} {hsa_resnames[idx]}: {rmsf_hsa[idx]:.2f} Å")

# --- 6. Linker Rg and Lys20-LNK distance ---
print("\n=== Linker Dynamics ===")

u3 = mda.Universe(PRMTOP, DCD)
lnk_sel = u3.select_atoms('resid 129')
lys_nz2 = u3.select_atoms('resid 117 and name NZ')
lnk_c2 = u3.select_atoms('resid 129 and name C')

lnk_rg = []
lys_lnk_c = []

for ts in u3.trajectory[::STEP]:
    lnk_rg.append(lnk_sel.radius_of_gyration())
    lys_lnk_c.append(np.linalg.norm(lys_nz2.positions[0] - lnk_c2.positions[0]))

lnk_rg = np.array(lnk_rg)
lys_lnk_c = np.array(lys_lnk_c)

print(f"  LNK Rg: mean={np.mean(lnk_rg):.2f} ± {np.std(lnk_rg):.2f} Å")
print(f"  Lys20 NZ - LNK C: mean={np.mean(lys_lnk_c):.3f} ± {np.std(lys_lnk_c):.3f} Å")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax1.plot(time_ns, lys_lnk_c, 'b-', alpha=0.7, label=f'Lys20 NZ-LNK C ({np.mean(lys_lnk_c):.3f} Å)', linewidth=0.8)
ax1.axhline(y=1.5, color='k', linestyle='--', alpha=0.3, label='Amide bond ref (1.5 Å)')
ax1.set_ylabel('Bond Length (Å)')
ax1.set_title('ECD v2 Lys20-LNK Amide Bond Length')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.plot(time_ns, lnk_rg, 'purple', alpha=0.7, label=f'LNK Rg ({np.mean(lnk_rg):.2f} Å)', linewidth=0.8)
ax2.set_xlabel('Time (ns)')
ax2.set_ylabel('Radius of Gyration (Å)')
ax2.set_title('ECD v2 Linker Radius of Gyration')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/linker_dynamics.png", dpi=150)
print(f"  Saved: {OUTDIR}/linker_dynamics.png")

# --- 7. Thermodynamics ---
print("\n=== Thermodynamics (from prod.log) ===")
log_path = f"{EXP_F}/md/ecd_v2/prod.log"
if os.path.exists(log_path):
    prod = pd.read_csv(log_path, skiprows=1, header=None)
    prod.columns = ['Step', 'Time_ps', 'PE', 'KE', 'TE', 'Temp', 'Vol', 'Dens', 'Speed']
    print(f"  Frames: {len(prod)}")
    print(f"  Temp: {prod['Temp'].mean():.2f} ± {prod['Temp'].std():.2f} K")
    print(f"  Density: {prod['Dens'].mean():.4f} ± {prod['Dens'].std():.4f} g/mL")
    print(f"  Volume: {prod['Vol'].mean():.2f} ± {prod['Vol'].std():.2f} nm³")
    
    pe_first = np.mean(prod['PE'][:100])
    pe_last = np.mean(prod['PE'][-100:])
    print(f"  PE drift (first 100 vs last 100): {pe_last - pe_first:.1f} kJ/mol")
    
    pe_diff = np.abs(np.diff(prod['PE'].values))
    big_jumps = np.sum(pe_diff > 1000)
    print(f"  Large PE jumps (>1000 kJ/mol): {big_jumps} ({big_jumps/len(prod)*100:.1f}% of frames)")

# --- 8. Summary ---
print("\n=== Summary ===")
summary = {
    'total_frames': n_frames,
    'analyzed_frames': n_analyzed,
    'total_ns': n_frames * 10 / 1000,
    'bond_lys20_lnk_c_mean_A': float(np.mean(lys_lnk_c)),
    'bond_lys20_lnk_c_std_A': float(np.std(lys_lnk_c)),
    'ecd_internal_rmsd_mean_A': float(np.mean(rmsd_ecd_int)),
    'ecd_internal_rmsd_max_A': float(np.max(rmsd_ecd_int)),
    'peptide_internal_rmsd_mean_A': float(np.mean(rmsd_pep_int)),
    'peptide_internal_rmsd_max_A': float(np.max(rmsd_pep_int)),
    'linker_internal_rmsd_mean_A': float(np.mean(rmsd_lnk_int)),
    'linker_internal_rmsd_max_A': float(np.max(rmsd_lnk_int)),
    'hsa_internal_rmsd_mean_A': float(np.mean(rmsd_hsa_int)),
    'hsa_internal_rmsd_max_A': float(np.max(rmsd_hsa_int)),
    'ecd_hsa_com_mean_A': float(np.mean(ecd_hsa_dist)),
    'ecd_hsa_com_std_A': float(np.std(ecd_hsa_dist)),
    'pep_hsa_com_mean_A': float(np.mean(pep_hsa_dist)),
    'pep_hsa_com_std_A': float(np.std(pep_hsa_dist)),
    'lnk_hsa_com_mean_A': float(np.mean(lnk_hsa_dist)),
    'lnk_hsa_com_std_A': float(np.std(lnk_hsa_dist)),
    'ecd_pep_com_mean_A': float(np.mean(ecd_pep_dist)),
    'ecd_pep_com_std_A': float(np.std(ecd_pep_dist)),
    'pep_lnk_com_mean_A': float(np.mean(pep_lnk_dist)),
    'pep_lnk_com_std_A': float(np.std(pep_lnk_dist)),
    'ecd_rmsf_mean_A': float(np.mean(rmsf_ecd)),
    'ecd_rmsf_max_A': float(np.max(rmsf_ecd)),
    'peptide_rmsf_mean_A': float(np.mean(rmsf_pep)),
    'peptide_rmsf_max_A': float(np.max(rmsf_pep)),
    'hsa_rmsf_mean_A': float(np.mean(rmsf_hsa)),
    'hsa_rmsf_max_A': float(np.max(rmsf_hsa)),
    'lnk_rg_mean_A': float(np.mean(lnk_rg)),
    'lnk_rg_std_A': float(np.std(lnk_rg)),
}

pd.DataFrame([summary]).to_csv(f"{OUTDIR}/summary_corrected.csv", index=False)
print(f"  Saved: {OUTDIR}/summary_corrected.csv")

for k, v in summary.items():
    if isinstance(v, float):
        print(f"  {k:40s}: {v:.3f}")
    else:
        print(f"  {k:40s}: {v}")

print(f"\nAnalysis complete. Output: {OUTDIR}/")
