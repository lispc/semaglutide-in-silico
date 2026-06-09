#!/usr/bin/env python3
"""
Full 100 ns analysis of ECD v2 production MD.

Corrected residue numbering based on system_ecd_v2.prmtop:
- ECD (GLP-1R): resid 1-100
- Peptide (semaglutide): resid 101-128
- Linker (LNK): resid 129
- HSA: resid 130-711

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

# --- 1. Internal RMSD (self-aligned) ---
print("\n=== Internal RMSD (self-aligned) ===")

# Load once with in_memory for speed
u_main = mda.Universe(PRMTOP, DCD)
ref_main = u_main.copy()
ref_main.trajectory[0]

# Align all protein heavy atoms once
print("  Aligning trajectory (in memory)...")
align.AlignTraj(u_main, ref_main, select='protein and name CA', in_memory=True).run(step=STEP)

# Now compute RMSD for each component using the aligned trajectory
sel_ecd = u_main.select_atoms('name CA and resid 1-100')
sel_pep = u_main.select_atoms('name CA and resid 101-128')
sel_lnk = u_main.select_atoms('resid 129 and not name H*')
sel_hsa = u_main.select_atoms('name CA and resid 130-711')

def compute_rmsd(sel, ref_sel):
    vals = []
    for ts in u_main.trajectory[::STEP]:
        vals.append(rms.rmsd(sel.positions, ref_sel.positions, superposition=False))
    return np.array(vals)

rmsd_ecd = compute_rmsd(sel_ecd, ref_main.select_atoms('name CA and resid 1-100'))
rmsd_pep = compute_rmsd(sel_pep, ref_main.select_atoms('name CA and resid 101-128'))
rmsd_lnk = compute_rmsd(sel_lnk, ref_main.select_atoms('resid 129 and not name H*'))
rmsd_hsa = compute_rmsd(sel_hsa, ref_main.select_atoms('name CA and resid 130-711'))

for label, vals in [("ECD", rmsd_ecd), ("Peptide", rmsd_pep), ("Linker", rmsd_lnk), ("HSA", rmsd_hsa)]:
    n = len(vals)
    drift = np.mean(vals[-n//10:]) - np.mean(vals[:n//10])
    print(f"  {label:12s}: mean={np.mean(vals):.3f} Å  max={np.max(vals):.3f} Å  drift={drift:+.3f} Å")

time_ns = np.arange(len(rmsd_ecd)) * (n_frames * 10 / 1000) / len(rmsd_ecd)

fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(time_ns, rmsd_ecd, 'b-', alpha=0.7, label=f'ECD ({np.mean(rmsd_ecd):.2f} Å)', linewidth=0.8)
ax.plot(time_ns, rmsd_pep, 'g-', alpha=0.7, label=f'Peptide ({np.mean(rmsd_pep):.2f} Å)', linewidth=0.8)
ax.plot(time_ns, rmsd_lnk, 'purple', alpha=0.7, label=f'Linker ({np.mean(rmsd_lnk):.2f} Å)', linewidth=0.8)
ax.plot(time_ns, rmsd_hsa, 'r-', alpha=0.7, label=f'HSA ({np.mean(rmsd_hsa):.2f} Å)', linewidth=0.8)
ax.set_xlabel('Time (ns)')
ax.set_ylabel('RMSD (Å)')
ax.set_title('ECD v2 Internal RMSD (100 ns, self-aligned)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsd_internal.png", dpi=150)
print(f"  Saved: {OUTDIR}/rmsd_internal.png")

# Save RMSD data
pd.DataFrame({
    'time_ns': time_ns,
    'ecd': rmsd_ecd,
    'peptide': rmsd_pep,
    'linker': rmsd_lnk,
    'hsa': rmsd_hsa,
}).to_csv(f"{OUTDIR}/rmsd_internal.csv", index=False)

# --- 2. RMSD Plateau Check ---
print("\n=== RMSD Plateau Analysis ===")

def check_plateau(vals, window=50):
    n = len(vals)
    if n < 2*window:
        window = n // 4
    early = np.mean(vals[:window])
    mid = np.mean(vals[n//2-window//2:n//2+window//2])
    late = np.mean(vals[-window:])
    return early, mid, late

for label, vals in [("ECD", rmsd_ecd), ("Peptide", rmsd_pep), ("Linker", rmsd_lnk), ("HSA", rmsd_hsa)]:
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

# --- 3. Internal RMSF ---
print("\n=== Internal RMSF (self-aligned) ===")

# Use the already-aligned trajectory from RMSD step
ecd_ca = u_main.select_atoms('name CA and resid 1-100')
pep_ca = u_main.select_atoms('name CA and resid 101-128')
hsa_ca = u_main.select_atoms('name CA and resid 130-711')
lnk_heavy = u_main.select_atoms('resid 129 and not name H*')

r_ecd = RMSF(ecd_ca).run(step=STEP)
rmsf_ecd = r_ecd.results.rmsf
ecd_resids = ecd_ca.resids
ecd_resnames = ecd_ca.resnames
print(f"  {'ECD':12s}: mean={np.mean(rmsf_ecd):.3f} Å  max={np.max(rmsf_ecd):.3f} Å")

r_pep = RMSF(pep_ca).run(step=STEP)
rmsf_pep = r_pep.results.rmsf
pep_resids = pep_ca.resids
pep_resnames = pep_ca.resnames
print(f"  {'Peptide':12s}: mean={np.mean(rmsf_pep):.3f} Å  max={np.max(rmsf_pep):.3f} Å")

r_hsa = RMSF(hsa_ca).run(step=STEP)
rmsf_hsa = r_hsa.results.rmsf
hsa_resids = hsa_ca.resids
hsa_resnames = hsa_ca.resnames
print(f"  {'HSA':12s}: mean={np.mean(rmsf_hsa):.3f} Å  max={np.max(rmsf_hsa):.3f} Å")

r_lnk = RMSF(lnk_heavy).run(step=STEP)
rmsf_lnk = r_lnk.results.rmsf
print(f"  {'Linker':12s}: mean={np.mean(rmsf_lnk):.3f} Å  max={np.max(rmsf_lnk):.3f} Å")

fig, ax = plt.subplots(figsize=(14, 5))
ax.plot(ecd_resids, rmsf_ecd, 'bo-', alpha=0.5, markersize=2, linewidth=0.5, label='ECD')
ax.plot(pep_resids, rmsf_pep, 'go-', alpha=0.7, markersize=4, linewidth=1, label='Peptide')
ax.plot(hsa_resids, rmsf_hsa, 'ro-', alpha=0.5, markersize=2, linewidth=0.5, label='HSA')
ax.axvline(x=100.5, color='gray', linestyle='--', alpha=0.5)
ax.axvline(x=128.5, color='gray', linestyle='--', alpha=0.5)
ax.axvline(x=129.5, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel('Residue ID')
ax.set_ylabel('RMSF (Å)')
ax.set_title('ECD v2 RMSF per Residue (100 ns, internal)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsf_internal.png", dpi=150)
print(f"  Saved: {OUTDIR}/rmsf_internal.png")

# --- 4. HSA High RMSF Residues ---
print("\n=== HSA High RMSF Residues (>5 Å) ===")
for i, (rid, rnm, val) in enumerate(zip(hsa_resids, hsa_resnames, rmsf_hsa)):
    if val > 5.0:
        print(f"  Residue {rid} {rnm}: {val:.2f} Å")

# Top 10
print("\n=== HSA Top 10 Highest RMSF ===")
top10 = np.argsort(rmsf_hsa)[-10:][::-1]
for idx in top10:
    print(f"  Residue {hsa_resids[idx]} {hsa_resnames[idx]}: {rmsf_hsa[idx]:.2f} Å")

# --- 5. Peptide-Linker Interface ---
print("\n=== Peptide-Linker Interface ===")

u3 = mda.Universe(PRMTOP, DCD)
ref3 = u3.copy()
ref3.trajectory[0]

# Key distances
lys_nz = u3.select_atoms('resid 117 and name NZ')  # Lys20 NZ in semaglutide
lnk_com = u3.select_atoms('resid 129')              # Linker center of mass

lys_lnk_dist = []
lnk_rg = []

for ts in u3.trajectory[::STEP]:
    d1 = np.linalg.norm(lys_nz.positions[0] - lnk_com.center_of_mass())
    lys_lnk_dist.append(d1)
    rg = lnk_com.radius_of_gyration()
    lnk_rg.append(rg)

lys_lnk_dist = np.array(lys_lnk_dist)
lnk_rg = np.array(lnk_rg)

print(f"  Lys20 NZ - LNK COM distance: {np.mean(lys_lnk_dist):.3f} ± {np.std(lys_lnk_dist):.3f} Å")
print(f"  LNK radius of gyration:      {np.mean(lnk_rg):.2f} ± {np.std(lnk_rg):.2f} Å")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax1.plot(time_ns, lys_lnk_dist, 'b-', alpha=0.7, label=f'Lys20 NZ-LNK COM ({np.mean(lys_lnk_dist):.2f} Å)', linewidth=0.8)
ax1.axhline(y=15.0, color='k', linestyle='--', alpha=0.3, label='Close contact ref (15 Å)')
ax1.set_ylabel('Distance (Å)')
ax1.set_title('ECD v2 Lys20 - Linker Distance')
ax1.legend()
ax1.grid(True, alpha=0.3)

ax2.plot(time_ns, lnk_rg, 'purple', alpha=0.7, label=f'LNK Rg ({np.mean(lnk_rg):.2f} Å)', linewidth=0.8)
ax2.set_xlabel('Time (ns)')
ax2.set_ylabel('Radius of Gyration (Å)')
ax2.set_title('ECD v2 Linker Radius of Gyration')
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/peptide_linker_interface.png", dpi=150)
print(f"  Saved: {OUTDIR}/peptide_linker_interface.png")

# --- 6. Thermodynamics ---
print("\n=== Thermodynamics (from prod.log) ===")
log_path = f"{EXP_F}/md/ecd_v2/prod.log"
if os.path.exists(log_path):
    prod = pd.read_csv(log_path, skiprows=1, header=None)
    prod.columns = ['Step', 'Time_ps', 'PE', 'KE', 'TE', 'Temp', 'Vol', 'Dens', 'Speed']
    print(f"  Frames: {len(prod)}")
    print(f"  Temp: {prod['Temp'].mean():.2f} ± {prod['Temp'].std():.2f} K")
    print(f"  Density: {prod['Dens'].mean():.4f} ± {prod['Dens'].std():.4f} g/mL")
    print(f"  Volume: {prod['Vol'].mean():.2f} ± {prod['Vol'].std():.2f} nm³")
    
    # PE drift
    pe_first = np.mean(prod['PE'][:100])
    pe_last = np.mean(prod['PE'][-100:])
    print(f"  PE drift (first 100 vs last 100): {pe_last - pe_first:.1f} kJ/mol")
    
    # Large jumps
    pe_diff = np.abs(np.diff(prod['PE'].values))
    big_jumps = np.sum(pe_diff > 1000)
    print(f"  Large PE jumps (>1000 kJ/mol): {big_jumps} ({big_jumps/len(prod)*100:.1f}% of frames)")

# --- 7. Summary ---
print("\n=== Summary ===")
summary = {
    'total_frames': n_frames,
    'analyzed_frames': n_analyzed,
    'total_ns': n_frames * 10 / 1000,
    'ecd_rmsd_mean_A': float(np.mean(rmsd_ecd)),
    'ecd_rmsd_max_A': float(np.max(rmsd_ecd)),
    'ecd_rmsd_drift_A': float(np.mean(rmsd_ecd[-len(rmsd_ecd)//10:]) - np.mean(rmsd_ecd[:len(rmsd_ecd)//10])),
    'peptide_rmsd_mean_A': float(np.mean(rmsd_pep)),
    'peptide_rmsd_max_A': float(np.max(rmsd_pep)),
    'peptide_rmsd_drift_A': float(np.mean(rmsd_pep[-len(rmsd_pep)//10:]) - np.mean(rmsd_pep[:len(rmsd_pep)//10])),
    'linker_rmsd_mean_A': float(np.mean(rmsd_lnk)),
    'linker_rmsd_max_A': float(np.max(rmsd_lnk)),
    'linker_rmsd_drift_A': float(np.mean(rmsd_lnk[-len(rmsd_lnk)//10:]) - np.mean(rmsd_lnk[:len(rmsd_lnk)//10])),
    'hsa_rmsd_mean_A': float(np.mean(rmsd_hsa)),
    'hsa_rmsd_max_A': float(np.max(rmsd_hsa)),
    'hsa_rmsd_drift_A': float(np.mean(rmsd_hsa[-len(rmsd_hsa)//10:]) - np.mean(rmsd_hsa[:len(rmsd_hsa)//10])),
    'ecd_rmsf_mean_A': float(np.mean(rmsf_ecd)),
    'ecd_rmsf_max_A': float(np.max(rmsf_ecd)),
    'peptide_rmsf_mean_A': float(np.mean(rmsf_pep)),
    'peptide_rmsf_max_A': float(np.max(rmsf_pep)),
    'hsa_rmsf_mean_A': float(np.mean(rmsf_hsa)),
    'hsa_rmsf_max_A': float(np.max(rmsf_hsa)),
    'lys20_lnk_com_mean_A': float(np.mean(lys_lnk_dist)),
    'lys20_lnk_com_std_A': float(np.std(lys_lnk_dist)),
    'lnk_rg_mean_A': float(np.mean(lnk_rg)),
    'lnk_rg_std_A': float(np.std(lnk_rg)),
}

pd.DataFrame([summary]).to_csv(f"{OUTDIR}/summary.csv", index=False)
print(f"  Saved: {OUTDIR}/summary.csv")

for k, v in summary.items():
    if isinstance(v, float):
        print(f"  {k:35s}: {v:.3f}")
    else:
        print(f"  {k:35s}: {v}")

print(f"\nAnalysis complete. Output: {OUTDIR}/")
