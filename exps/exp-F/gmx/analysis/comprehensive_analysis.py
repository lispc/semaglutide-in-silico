#!/usr/bin/env python3
"""
Comprehensive GROMACS membrane MD analysis script.
Designed for 200 ns production run validation.

Usage:
    cd exps/exp-F/gmx && python analysis/comprehensive_analysis.py

Selections (GROMACS gro/xtc, resindex-based due to multi-chain renumbering):
  ECD:      resindex 0-1112   (1113 residues, 17709 atoms)
  GLP-1:    resindex 1113-1138 (26 residues, 390 atoms)
  LNK:      resindex 1139     (1 residue, 47 atoms)
  Lipid:    resname PA PC OL CHL
  Water:    resname WAT
"""

import os
import sys
import subprocess
import json
import warnings
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.pca import PCA as MDA_PCA
from MDAnalysis.analysis.rdf import InterRDF

warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

# =============================================================================
# CONFIGURATION
# =============================================================================
REPO = "/home/scroll/personal/semaglutide-in-silico"
GMX_DIR = f"{REPO}/exps/exp-F/gmx"
OUTDIR = f"{GMX_DIR}/analysis"
GRO = f"{GMX_DIR}/npt.gro"
XTC = f"{GMX_DIR}/md.xtc"
EDR = f"{GMX_DIR}/md.edr"

# Trajectory sampling step (every Nth frame). XTC saved every 2 ps, so step=10 = 20 ps.
STEP = 10

os.makedirs(OUTDIR, exist_ok=True)

# Energy terms to extract from EDR
ENERGY_TERMS = [
    "Potential", "Kinetic-En.", "Total-Energy", "Temperature",
    "Pressure", "Volume", "Density",
    "Coulomb-14", "LJ-14", "Coulomb-(SR)", "LJ-(SR)",
    "Coul.-recip.", "LJ-recip.", "Constraint-Viol.",
    "Box-X", "Box-Y", "Box-Z",
]

# =============================================================================
# LOGGING
# =============================================================================
log_file = open(f"{OUTDIR}/analysis.log", "w")

def log(msg):
    print(msg)
    log_file.write(msg + "\n")
    log_file.flush()

log("="*70)
log("COMPREHENSIVE GROMACS MEMBRANE MD ANALYSIS")
log("="*70)
log(f"Working directory: {GMX_DIR}")
log(f"Output directory:  {OUTDIR}")
log(f"Trajectory:        {XTC}")
log(f"Topology:          {GRO}")
log(f"Energy file:       {EDR}")
log(f"Analysis step:     every {STEP} frames")

# =============================================================================
# 1. LOAD UNIVERSE & DEFINE SELECTIONS
# =============================================================================
log("\n[1] Loading universe...")
u = mda.Universe(GRO, XTC)
n_frames = len(u.trajectory)
n_analyzed = (n_frames + STEP - 1) // STEP
log(f"  Total frames: {n_frames}")
log(f"  Analyzing:    {n_analyzed} (step={STEP})")
log(f"  Total atoms:  {len(u.atoms)}")

# Create reference universe for alignment
ref = mda.Universe(GRO)

# Selections
sel_protein = u.select_atoms('protein')                    # ECD + GLP-1 (excludes LNK)
sel_ecd     = u.select_atoms('resindex 0-1112')            # ECD only
sel_glp1    = u.select_atoms('resindex 1113-1138')         # GLP-1 only
sel_lnk     = u.select_atoms('resindex 1139')              # Linker only
sel_lipid   = u.select_atoms('resname PA PC OL CHL')
sel_water   = u.select_atoms('resname WAT')
sel_ions    = u.select_atoms('resname K+ Cl-')

# CA selections for backbone analysis
sel_protein_ca = sel_protein.select_atoms('name CA')
sel_ecd_ca     = sel_ecd.select_atoms('name CA')
sel_glp1_ca    = sel_glp1.select_atoms('name CA')

log(f"  Protein (protein):     {len(sel_protein.atoms):6d} atoms, {len(sel_protein.residues):4d} res")
log(f"  ECD (resindex 0-1112): {len(sel_ecd.atoms):6d} atoms, {len(sel_ecd.residues):4d} res")
log(f"  GLP-1 (resindex 1113-1138): {len(sel_glp1.atoms):6d} atoms, {len(sel_glp1.residues):4d} res")
log(f"  LNK (resindex 1139):   {len(sel_lnk.atoms):6d} atoms, {len(sel_lnk.residues):4d} res")
log(f"  Lipid:                 {len(sel_lipid.atoms):6d} atoms, {len(sel_lipid.residues):4d} res")
log(f"  Water:                 {len(sel_water.atoms):6d} atoms, {len(sel_water.residues):4d} res")
log(f"  Ions:                  {len(sel_ions.atoms):6d} atoms, {len(sel_ions.residues):4d} res")

# Ref selections
ref_protein_ca = ref.select_atoms('protein and name CA')
ref_ecd_ca     = ref.select_atoms('resindex 0-1112 and name CA')
ref_glp1_ca    = ref.select_atoms('resindex 1113-1138 and name CA')
ref_lnk        = ref.select_atoms('resindex 1139')

# =============================================================================
# HELPER: PBC unwrap for selected group
# =============================================================================
def unwrap_group(ag):
    """Unwrap atom group across PBC if possible."""
    if len(ag) > 0 and hasattr(ag, 'unwrap'):
        try:
            ag.unwrap()
        except Exception:
            pass

# =============================================================================
# 2. ENERGY ANALYSIS (via gmx energy)
# =============================================================================
log("\n[2] Extracting energy terms from EDR...")

def extract_energy(edr_path, terms, out_csv):
    """Use gmx energy to extract terms."""
    # Build xvg command
    terms_str = ' '.join(terms)
    cmd = f"echo '{terms_str}' | gmx energy -f {edr_path} -o {out_csv.replace('.csv','.xvg')} 2>&1"
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=60)
    except Exception as e:
        log(f"  Warning: gmx energy failed: {e}")
        return None
    
    # Parse xvg to csv
    if not os.path.exists(out_csv.replace('.csv','.xvg')):
        log(f"  Warning: xvg not created")
        return None
    
    data = {}
    with open(out_csv.replace('.csv','.xvg')) as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            if line.strip():
                parts = line.split()
                if len(parts) >= 2:
                    t = float(parts[0])
                    if 'time' not in data:
                        data['time'] = []
                        for _ in range(len(parts)-1):
                            data[f'col_{len(data)}'] = []
                    data['time'].append(t)
                    for i, val in enumerate(parts[1:]):
                        key = list(data.keys())[i+1]
                        data[key].append(float(val))
    
    if data:
        df = pd.DataFrame(data)
        df.to_csv(out_csv, index=False)
        log(f"  Saved: {out_csv}")
        return df
    return None

# Alternative: use panedr if available
try:
    import panedr
    log("  Using panedr for fast EDR parsing...")
    df_energy = panedr.edr_to_df(EDR)
    # Save subset
    avail = [c for c in ENERGY_TERMS if c in df_energy.columns]
    df_energy[['Time'] + avail].to_csv(f"{OUTDIR}/energy.csv", index=False)
    log(f"  Saved: {OUTDIR}/energy.csv (cols: {avail})")
    energy_ok = True
except ImportError:
    log("  panedr not available, skipping energy CSV (will use xvg directly)")
    energy_ok = False
    df_energy = None

# Quick energy plot if available
if energy_ok and 'Potential' in df_energy.columns:
    fig, axes = plt.subplots(3, 2, figsize=(12, 10))
    axes = axes.flatten()
    
    time_ns = df_energy['Time'].values / 1000.0
    plots = [
        ('Potential', 'Potential Energy (kJ/mol)'),
        ('Kinetic-En.', 'Kinetic Energy (kJ/mol)'),
        ('Total-Energy', 'Total Energy (kJ/mol)'),
        ('Temperature', 'Temperature (K)'),
        ('Pressure', 'Pressure (bar)'),
        ('Density', 'Density (kg/m³)'),
    ]
    for ax, (col, label) in zip(axes, plots):
        if col in df_energy.columns:
            ax.plot(time_ns, df_energy[col].values, linewidth=0.5)
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel(label)
            ax.set_title(label)
            ax.axvline(x=time_ns[-1], color='r', linestyle='--', alpha=0.3, label=f'End: {time_ns[-1]:.1f} ns')
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/energy_overview.png", dpi=150)
    plt.close()
    log(f"  Saved: {OUTDIR}/energy_overview.png")

# =============================================================================
# 3. RMSD ANALYSIS
# =============================================================================
log("\n[3] RMSD analysis...")

def compute_rmsd(u_obj, ref_obj, select_str, ref_select_str, label, step=STEP):
    """Compute RMSD for a selection."""
    log(f"  Computing RMSD: {label}...")
    R = rms.RMSD(u_obj, ref_obj,
                 select=select_str,
                 groupselections=None)
    R.run(step=step)
    data = R.results.rmsd
    df = pd.DataFrame(data, columns=['Frame', 'Time (ps)', 'RMSD'])
    df['Time (ns)'] = df['Time (ps)'] / 1000.0
    df = df.rename(columns={'RMSD': label})
    return df

# Align on ECD CA, then compute RMSDs for all components
# First: overall protein CA RMSD (no alignment, just raw)
rmsd_protein_raw = compute_rmsd(u, ref, 'protein and name CA', 'protein and name CA', 'Protein_CA_Raw', step=STEP)

# Align on ECD CA and recompute
log("  Aligning trajectory on ECD CA...")
# We need to iterate manually for aligned RMSD
rmsd_ecd_aligned = []
rmsd_glp1_aligned = []
rmsd_lnk_aligned = []
rmsd_protein_aligned = []
times = []

ref.trajectory[0]
ref_ecd_pos = ref_ecd_ca.positions.copy()
ref_glp1_pos = ref_glp1_ca.positions.copy()
ref_lnk_pos = ref_lnk.positions.copy()
ref_prot_pos = ref_protein_ca.positions.copy()

for ts in u.trajectory[::STEP]:
    # Unwrap ECD
    unwrap_group(sel_ecd)
    
    # Align on ECD CA
    align.alignto(sel_ecd_ca, ref_ecd_ca, select='name CA', weights='mass')
    
    # Compute RMSDs after alignment
    rmsd_ecd_aligned.append(rms.rmsd(sel_ecd_ca.positions, ref_ecd_pos, superposition=False))
    rmsd_glp1_aligned.append(rms.rmsd(sel_glp1_ca.positions, ref_glp1_pos, superposition=False))
    rmsd_lnk_aligned.append(rms.rmsd(sel_lnk.positions, ref_lnk_pos, superposition=False))
    rmsd_protein_aligned.append(rms.rmsd(sel_protein_ca.positions, ref_prot_pos, superposition=False))
    times.append(ts.time / 1000.0)

times = np.array(times)
rmsd_ecd_aligned = np.array(rmsd_ecd_aligned)
rmsd_glp1_aligned = np.array(rmsd_glp1_aligned)
rmsd_lnk_aligned = np.array(rmsd_lnk_aligned)
rmsd_protein_aligned = np.array(rmsd_protein_aligned)

df_rmsd = pd.DataFrame({
    'Time_ns': times,
    'RMSD_ECD_CA_A': rmsd_ecd_aligned * 10,  # nm -> Å
    'RMSD_GLP1_CA_A': rmsd_glp1_aligned * 10,
    'RMSD_LNK_Heavy_A': rmsd_lnk_aligned * 10,
    'RMSD_Protein_CA_A': rmsd_protein_aligned * 10,
})
df_rmsd.to_csv(f"{OUTDIR}/rmsd.csv", index=False)
log(f"  Saved: {OUTDIR}/rmsd.csv")

# Plot RMSD
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes[0,0].plot(times, df_rmsd['RMSD_Protein_CA_A'], 'k-', linewidth=0.5)
axes[0,0].set_title('Protein CA RMSD (Å)')
axes[0,0].set_ylabel('RMSD (Å)')
axes[0,1].plot(times, df_rmsd['RMSD_ECD_CA_A'], 'b-', linewidth=0.5)
axes[0,1].set_title('ECD CA RMSD (Å)')
axes[0,1].set_ylabel('RMSD (Å)')
axes[1,0].plot(times, df_rmsd['RMSD_GLP1_CA_A'], 'g-', linewidth=0.5)
axes[1,0].set_title('GLP-1 CA RMSD (Å)')
axes[1,0].set_xlabel('Time (ns)')
axes[1,0].set_ylabel('RMSD (Å)')
axes[1,1].plot(times, df_rmsd['RMSD_LNK_Heavy_A'], 'purple', linewidth=0.5)
axes[1,1].set_title('LNK Heavy Atom RMSD (Å)')
axes[1,1].set_xlabel('Time (ns)')
axes[1,1].set_ylabel('RMSD (Å)')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsd_overview.png", dpi=150)
plt.close()
log(f"  Saved: {OUTDIR}/rmsd_overview.png")

# Summary statistics
log(f"  RMSD stats (Å):")
log(f"    Protein CA:  {df_rmsd['RMSD_Protein_CA_A'].mean():.2f} ± {df_rmsd['RMSD_Protein_CA_A'].std():.2f}")
log(f"    ECD CA:      {df_rmsd['RMSD_ECD_CA_A'].mean():.2f} ± {df_rmsd['RMSD_ECD_CA_A'].std():.2f}")
log(f"    GLP-1 CA:    {df_rmsd['RMSD_GLP1_CA_A'].mean():.2f} ± {df_rmsd['RMSD_GLP1_CA_A'].std():.2f}")
log(f"    LNK Heavy:   {df_rmsd['RMSD_LNK_Heavy_A'].mean():.2f} ± {df_rmsd['RMSD_LNK_Heavy_A'].std():.2f}")

# =============================================================================
# 4. RADIUS OF GYRATION
# =============================================================================
log("\n[4] Radius of Gyration...")

def compute_rg(u_obj, selection, step=STEP):
    """Compute Rg for a selection."""
    rg_vals = []
    times = []
    for ts in u_obj.trajectory[::step]:
        unwrap_group(selection)
        rg_vals.append(selection.radius_of_gyration())
        times.append(ts.time / 1000.0)
    return np.array(times), np.array(rg_vals)

t_rg, rg_protein = compute_rg(u, sel_protein, STEP)
_, rg_ecd = compute_rg(u, sel_ecd, STEP)
_, rg_glp1 = compute_rg(u, sel_glp1, STEP)

df_rg = pd.DataFrame({
    'Time_ns': t_rg,
    'Rg_Protein_nm': rg_protein,
    'Rg_ECD_nm': rg_ecd,
    'Rg_GLP1_nm': rg_glp1,
    'Rg_Protein_A': rg_protein * 10,
    'Rg_ECD_A': rg_ecd * 10,
    'Rg_GLP1_A': rg_glp1 * 10,
})
df_rg.to_csv(f"{OUTDIR}/rg.csv", index=False)
log(f"  Saved: {OUTDIR}/rg.csv")

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
axes[0].plot(t_rg, df_rg['Rg_Protein_A'], 'k-', linewidth=0.5)
axes[0].set_title(f'Protein Rg: {df_rg["Rg_Protein_A"].mean():.1f}±{df_rg["Rg_Protein_A"].std():.1f} Å')
axes[0].set_ylabel('Rg (Å)')
axes[0].set_xlabel('Time (ns)')
axes[1].plot(t_rg, df_rg['Rg_ECD_A'], 'b-', linewidth=0.5)
axes[1].set_title(f'ECD Rg: {df_rg["Rg_ECD_A"].mean():.1f}±{df_rg["Rg_ECD_A"].std():.1f} Å')
axes[1].set_xlabel('Time (ns)')
axes[2].plot(t_rg, df_rg['Rg_GLP1_A'], 'g-', linewidth=0.5)
axes[2].set_title(f'GLP-1 Rg: {df_rg["Rg_GLP1_A"].mean():.1f}±{df_rg["Rg_GLP1_A"].std():.1f} Å')
axes[2].set_xlabel('Time (ns)')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rg_overview.png", dpi=150)
plt.close()
log(f"  Saved: {OUTDIR}/rg_overview.png")

# =============================================================================
# 5. RMSF ANALYSIS
# =============================================================================
log("\n[5] RMSF analysis...")

def compute_rmsf(u_obj, selection, ref_selection, label, step=STEP):
    """Compute RMSF aligned on the selection's CA."""
    log(f"  Computing RMSF: {label}...")
    # Collect coordinates
    coords = []
    for ts in u_obj.trajectory[::step]:
        unwrap_group(selection)
        align.alignto(selection, ref_selection, select='name CA', weights='mass')
        coords.append(selection.positions.copy())
    coords = np.array(coords)
    mean_pos = coords.mean(axis=0)
    rmsf = np.sqrt(((coords - mean_pos)**2).sum(axis=2)).mean(axis=0)
    return rmsf

# RMSF for ECD CA
rmsf_ecd = compute_rmsf(u, sel_ecd_ca, ref_ecd_ca, 'ECD CA', STEP)
# RMSF for GLP-1 CA
rmsf_glp1 = compute_rmsf(u, sel_glp1_ca, ref_glp1_ca, 'GLP-1 CA', STEP)

df_rmsf_ecd = pd.DataFrame({
    'ResIndex': range(len(sel_ecd_ca)),
    'ResID': [r.resid for r in sel_ecd_ca.residues],
    'ResName': [r.resname for r in sel_ecd_ca.residues],
    'RMSF_A': rmsf_ecd * 10,
})
df_rmsf_ecd.to_csv(f"{OUTDIR}/rmsf_ecd.csv", index=False)

df_rmsf_glp1 = pd.DataFrame({
    'ResIndex': range(len(sel_glp1_ca)),
    'ResID': [r.resid for r in sel_glp1_ca.residues],
    'ResName': [r.resname for r in sel_glp1_ca.residues],
    'RMSF_A': rmsf_glp1 * 10,
})
df_rmsf_glp1.to_csv(f"{OUTDIR}/rmsf_glp1.csv", index=False)
log(f"  Saved: {OUTDIR}/rmsf_ecd.csv, rmsf_glp1.csv")

# Plot RMSF
fig, axes = plt.subplots(2, 1, figsize=(14, 8))
axes[0].plot(df_rmsf_ecd['ResIndex'], df_rmsf_ecd['RMSF_A'], 'b-', linewidth=0.5)
axes[0].set_title('ECD CA RMSF (Å)')
axes[0].set_ylabel('RMSF (Å)')
axes[0].axhline(y=df_rmsf_ecd['RMSF_A'].mean(), color='r', linestyle='--', alpha=0.5, label=f'Mean: {df_rmsf_ecd["RMSF_A"].mean():.1f} Å')
axes[0].legend()

axes[1].plot(df_rmsf_glp1['ResIndex'], df_rmsf_glp1['RMSF_A'], 'g-o', linewidth=1, markersize=3)
for i, (ri, rn, val) in enumerate(zip(df_rmsf_glp1['ResIndex'], df_rmsf_glp1['ResName'], df_rmsf_glp1['RMSF_A'])):
    if val > df_rmsf_glp1['RMSF_A'].mean() + 1.5 * df_rmsf_glp1['RMSF_A'].std():
        axes[1].annotate(f'{rn}{ri+1}', xy=(ri, val), xytext=(ri, val+0.5), fontsize=7, ha='center')
axes[1].set_title('GLP-1 CA RMSF (Å)')
axes[1].set_ylabel('RMSF (Å)')
axes[1].set_xlabel('Residue Index')
axes[1].axhline(y=df_rmsf_glp1['RMSF_A'].mean(), color='r', linestyle='--', alpha=0.5, label=f'Mean: {df_rmsf_glp1["RMSF_A"].mean():.1f} Å')
axes[1].legend()

plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsf_overview.png", dpi=150)
plt.close()
log(f"  Saved: {OUTDIR}/rmsf_overview.png")

# =============================================================================
# 6. LNK POSITION & MEMBRANE ANALYSIS
# =============================================================================
log("\n[6] LNK position & membrane analysis...")

# Phosphorus atoms for membrane headgroups
phosphorus = u.select_atoms('name P31 or name P')
if len(phosphorus) == 0:
    # Try common lipid phosphorus names
    phosphorus = u.select_atoms('name P* and (resname PA PC OL)')

log(f"  Phosphorus atoms: {len(phosphorus)}")

lnk_com_z = []
lnk_head_z = []
lnk_tail_z = []
mem_upper_z = []
mem_lower_z = []
mem_com_z = []
lnk_tilt = []
dist_lys_nz_to_lnk = []
times_pos = []

# LNK head = first few atoms (amide region), tail = last few atoms (fatty acid end)
lnk_head_atoms = sel_lnk.select_atoms('name C13 or name C14 or name C15 or name N16')
lnk_tail_atoms = sel_lnk.select_atoms('name C48 or name C49 or name C50 or name C51 or name C52 or name C53 or name C54 or name C55')

# Lys NZ (Lys26 in GLP-1, resindex 1129)
lys_nz = u.select_atoms('resindex 1129 and name NZ')

for ts in u.trajectory[::STEP]:
    # Unwrap protein and lipids
    unwrap_group(sel_protein)
    unwrap_group(sel_lipid)
    
    # Membrane phosphorus z-coords
    p_z = phosphorus.positions[:, 2]
    upper = np.percentile(p_z, 90)
    lower = np.percentile(p_z, 10)
    mem_com = p_z.mean()
    
    # LNK positions
    lnk_com = sel_lnk.center_of_mass()
    lnk_head_com = lnk_head_atoms.center_of_mass() if len(lnk_head_atoms) > 0 else lnk_com
    lnk_tail_com = lnk_tail_atoms.center_of_mass() if len(lnk_tail_atoms) > 0 else lnk_com
    
    # Tilt angle (head-to-tail vector vs z-axis)
    vec = lnk_tail_com - lnk_head_com
    vec_norm = np.linalg.norm(vec)
    tilt = np.degrees(np.arccos(abs(vec[2]) / vec_norm)) if vec_norm > 0 else 0
    
    # Lys NZ to LNK distance
    if len(lys_nz) > 0 and len(sel_lnk) > 0:
        d = np.linalg.norm(lys_nz.center_of_mass() - lnk_com)
        dist_lys_nz_to_lnk.append(d)
    else:
        dist_lys_nz_to_lnk.append(np.nan)
    
    lnk_com_z.append(lnk_com[2])
    lnk_head_z.append(lnk_head_com[2])
    lnk_tail_z.append(lnk_tail_com[2])
    mem_upper_z.append(upper)
    mem_lower_z.append(lower)
    mem_com_z.append(mem_com)
    lnk_tilt.append(tilt)
    times_pos.append(ts.time / 1000.0)

times_pos = np.array(times_pos)
lnk_com_z = np.array(lnk_com_z)
lnk_head_z = np.array(lnk_head_z)
lnk_tail_z = np.array(lnk_tail_z)
mem_upper_z = np.array(mem_upper_z)
mem_lower_z = np.array(mem_lower_z)
mem_com_z = np.array(mem_com_z)
lnk_tilt = np.array(lnk_tilt)
dist_lys_nz_to_lnk = np.array(dist_lys_nz_to_lnk)

df_pos = pd.DataFrame({
    'Time_ns': times_pos,
    'LNK_COM_z_nm': lnk_com_z,
    'LNK_Head_z_nm': lnk_head_z,
    'LNK_Tail_z_nm': lnk_tail_z,
    'Mem_Upper_z_nm': mem_upper_z,
    'Mem_Lower_z_nm': mem_lower_z,
    'Mem_COM_z_nm': mem_com_z,
    'LNK_Tilt_deg': lnk_tilt,
    'LysNZ_LNK_Dist_nm': dist_lys_nz_to_lnk,
    'LNK_COM_rel_mem_nm': lnk_com_z - mem_com_z,
    'LNK_Tail_rel_upper_nm': lnk_tail_z - mem_upper_z,
    'LNK_Head_rel_upper_nm': lnk_head_z - mem_upper_z,
})
df_pos.to_csv(f"{OUTDIR}/lnk_position.csv", index=False)
log(f"  Saved: {OUTDIR}/lnk_position.csv")

# Plot LNK position
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes[0,0].plot(times_pos, lnk_com_z * 10, 'purple', linewidth=0.5, label='LNK COM')
axes[0,0].plot(times_pos, mem_upper_z * 10, 'b--', linewidth=0.5, alpha=0.5, label='Membrane upper')
axes[0,0].plot(times_pos, mem_lower_z * 10, 'b--', linewidth=0.5, alpha=0.5, label='Membrane lower')
axes[0,0].set_title('LNK COM Z Position (Å)')
axes[0,0].set_ylabel('Z (Å)')
axes[0,0].legend(fontsize=7)

axes[0,1].plot(times_pos, (lnk_com_z - mem_com_z) * 10, 'purple', linewidth=0.5)
axes[0,1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[0,1].set_title('LNK COM relative to Membrane COM (Å)')
axes[0,1].set_ylabel('ΔZ (Å)')

axes[1,0].plot(times_pos, lnk_tilt, 'red', linewidth=0.5)
axes[1,0].set_title('LNK Tilt Angle (°)')
axes[1,0].set_ylabel('Tilt (°)')
axes[1,0].set_xlabel('Time (ns)')

if not np.all(np.isnan(dist_lys_nz_to_lnk)):
    axes[1,1].plot(times_pos, dist_lys_nz_to_lnk * 10, 'green', linewidth=0.5)
    axes[1,1].set_title('Lys NZ — LNK Distance (Å)')
    axes[1,1].set_ylabel('Distance (Å)')
    axes[1,1].set_xlabel('Time (ns)')
else:
    axes[1,1].text(0.5, 0.5, 'Lys NZ not found', ha='center', va='center', transform=axes[1,1].transAxes)

plt.tight_layout()
plt.savefig(f"{OUTDIR}/lnk_position.png", dpi=150)
plt.close()
log(f"  Saved: {OUTDIR}/lnk_position.png")

log(f"  LNK stats:")
log(f"    COM relative to mem: {(lnk_com_z - mem_com_z).mean()*10:.2f} ± {(lnk_com_z - mem_com_z).std()*10:.2f} Å")
log(f"    Tilt: {lnk_tilt.mean():.1f} ± {lnk_tilt.std():.1f}°")
if not np.all(np.isnan(dist_lys_nz_to_lnk)):
    log(f"    Lys-LNK dist: {np.nanmean(dist_lys_nz_to_lnk)*10:.2f} ± {np.nanstd(dist_lys_nz_to_lnk)*10:.2f} Å")

# =============================================================================
# 7. MEMBRANE PROPERTIES
# =============================================================================
log("\n[7] Membrane properties...")

# Area per lipid: box_xy / n_lipids_per_leaflet
# Membrane thickness: average distance between upper and lower leaflet P planes

membrane_thickness = (mem_upper_z - mem_lower_z) * 10  # nm -> Å
n_lipid_molecules = len(phosphorus.residues) if len(phosphorus) > 0 else len(sel_lipid.residues)

# Get box dimensions
df_box = None
if energy_ok and 'Box-X' in df_energy.columns:
    box_x = df_energy['Box-X'].values  # nm
    box_y = df_energy['Box-Y'].values if 'Box-Y' in df_energy.columns else box_x
    box_z = df_energy['Box-Z'].values if 'Box-Z' in df_energy.columns else None
    area_per_lipid = (box_x * box_y * 100) / (n_lipid_molecules / 2)  # nm² -> Å², divide by 2 leaflets
    
    df_mem = pd.DataFrame({
        'Time_ns': df_energy['Time'].values / 1000.0,
        'Box_X_nm': box_x,
        'Box_Y_nm': box_y,
        'Box_Z_nm': box_z if box_z is not None else np.nan,
        'Area_per_lipid_A2': area_per_lipid,
        'Membrane_thickness_A': np.interp(df_energy['Time'].values / 1000.0, times_pos, membrane_thickness),
    })
    df_mem.to_csv(f"{OUTDIR}/membrane_properties.csv", index=False)
    log(f"  Saved: {OUTDIR}/membrane_properties.csv")
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(df_mem['Time_ns'], df_mem['Area_per_lipid_A2'], 'b-', linewidth=0.5)
    axes[0].set_title(f'Area per Lipid: {df_mem["Area_per_lipid_A2"].mean():.1f}±{df_mem["Area_per_lipid_A2"].std():.1f} Å²')
    axes[0].set_ylabel('Area (Å²)')
    axes[0].set_xlabel('Time (ns)')
    
    axes[1].plot(times_pos, membrane_thickness, 'b-', linewidth=0.5)
    axes[1].set_title(f'Membrane Thickness: {membrane_thickness.mean():.1f}±{membrane_thickness.std():.1f} Å')
    axes[1].set_ylabel('Thickness (Å)')
    axes[1].set_xlabel('Time (ns)')
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/membrane_properties.png", dpi=150)
    plt.close()
    log(f"  Saved: {OUTDIR}/membrane_properties.png")
else:
    log("  Skipping membrane area (Box data unavailable)")
    df_mem = pd.DataFrame({'Time_ns': times_pos, 'Membrane_thickness_A': membrane_thickness})
    df_mem.to_csv(f"{OUTDIR}/membrane_thickness.csv", index=False)

# =============================================================================
# 8. HYDROGEN BOND ANALYSIS
# =============================================================================
log("\n[8] Hydrogen bond analysis...")

def run_hbonds(u_obj, sel1, sel2, label, step=STEP, d_a_cutoff=3.5, angle_cutoff=150):
    """Run hydrogen bond analysis between two selections."""
    log(f"  H-bonds: {label}...")
    try:
        hb = HydrogenBondAnalysis(
            universe=u_obj,
            donors_sel=sel1,
            hydrogens_sel=sel1,
            acceptors_sel=sel2,
            d_a_cutoff=d_a_cutoff,
            d_h_a_angle_cutoff=angle_cutoff,
            update_selections=False,
        )
        hb.run(step=step, verbose=False)
        
        counts = hb.results.hbonds.shape[0] if len(hb.results.hbonds) > 0 else 0
        # Per-frame counts
        frame_counts = {}
        if len(hb.results.hbonds) > 0:
            for row in hb.results.hbonds:
                frame = int(row[0])
                frame_counts[frame] = frame_counts.get(frame, 0) + 1
        
        # Align with our times array
        n_frames_analyzed = len(u_obj.trajectory[::step])
        counts_per_frame = [frame_counts.get(i, 0) for i in range(n_frames_analyzed)]
        return np.array(counts_per_frame)
    except Exception as e:
        log(f"    Warning: H-bond analysis failed: {e}")
        return None

# Protein internal H-bonds
hb_internal = run_hbonds(u, 'protein', 'protein', 'Protein internal', STEP)
# Protein-Lipid H-bonds
hb_prot_lipid = run_hbonds(u, 'protein', 'resname PA PC OL', 'Protein-Lipid', STEP)
# GLP-1-ECD H-bonds
hb_glp1_ecd = run_hbonds(u, 'resindex 1113-1138', 'resindex 0-1112', 'GLP-1-ECD', STEP)

if hb_internal is not None and hb_prot_lipid is not None:
    df_hb = pd.DataFrame({
        'Time_ns': times_pos[:len(hb_internal)],
        'HB_Internal': hb_internal,
        'HB_Protein_Lipid': hb_prot_lipid[:len(hb_internal)],
    })
    if hb_glp1_ecd is not None:
        df_hb['HB_GLP1_ECD'] = hb_glp1_ecd[:len(hb_internal)]
    df_hb.to_csv(f"{OUTDIR}/hbonds.csv", index=False)
    log(f"  Saved: {OUTDIR}/hbonds.csv")
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 6))
    axes[0].plot(df_hb['Time_ns'], df_hb['HB_Internal'], 'b-', linewidth=0.5)
    axes[0].set_title(f'Protein Internal H-Bonds: {df_hb["HB_Internal"].mean():.1f}±{df_hb["HB_Internal"].std():.1f}')
    axes[0].set_ylabel('Count')
    
    axes[1].plot(df_hb['Time_ns'], df_hb['HB_Protein_Lipid'], 'purple', linewidth=0.5)
    axes[1].set_title(f'Protein-Lipid H-Bonds: {df_hb["HB_Protein_Lipid"].mean():.1f}±{df_hb["HB_Protein_Lipid"].std():.1f}')
    axes[1].set_ylabel('Count')
    axes[1].set_xlabel('Time (ns)')
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/hbonds.png", dpi=150)
    plt.close()
    log(f"  Saved: {OUTDIR}/hbonds.png")
else:
    log("  H-bond analysis partially failed")

# =============================================================================
# 9. PROTEIN-LIPID CONTACTS
# =============================================================================
log("\n[9] Protein-lipid contacts...")

def compute_contacts(u_obj, sel1, sel2, cutoff_nm=0.5, step=STEP):
    """Compute number of contacts between two selections."""
    log(f"  Contacts: {sel1} - {sel2}...")
    counts = []
    times_c = []
    for ts in u_obj.trajectory[::step]:
        # Use minimum distance search
        pairs = mda.lib.distances.capped_distance(
            sel1.positions, sel2.positions,
            max_cutoff=cutoff_nm,
            return_distances=False,
            box=u_obj.dimensions,
        )
        counts.append(len(pairs))
        times_c.append(ts.time / 1000.0)
    return np.array(times_c), np.array(counts)

# ECD-Lipid contacts
t_c, cont_ecd_lipid = compute_contacts(u, sel_ecd, sel_lipid, 0.5, STEP)
# GLP-1-Lipid contacts
_, cont_glp1_lipid = compute_contacts(u, sel_glp1, sel_lipid, 0.5, STEP)
# LNK-Lipid contacts
_, cont_lnk_lipid = compute_contacts(u, sel_lnk, sel_lipid, 0.5, STEP)

df_cont = pd.DataFrame({
    'Time_ns': t_c,
    'Contacts_ECD_Lipid': cont_ecd_lipid,
    'Contacts_GLP1_Lipid': cont_glp1_lipid,
    'Contacts_LNK_Lipid': cont_lnk_lipid,
})
df_cont.to_csv(f"{OUTDIR}/contacts.csv", index=False)
log(f"  Saved: {OUTDIR}/contacts.csv")

fig, axes = plt.subplots(3, 1, figsize=(10, 8))
axes[0].plot(t_c, cont_ecd_lipid, 'b-', linewidth=0.5)
axes[0].set_title(f'ECD-Lipid Contacts (<5Å): {cont_ecd_lipid.mean():.1f}±{cont_ecd_lipid.std():.1f}')
axes[0].set_ylabel('Contacts')
axes[1].plot(t_c, cont_glp1_lipid, 'g-', linewidth=0.5)
axes[1].set_title(f'GLP-1-Lipid Contacts (<5Å): {cont_glp1_lipid.mean():.1f}±{cont_glp1_lipid.std():.1f}')
axes[1].set_ylabel('Contacts')
axes[2].plot(t_c, cont_lnk_lipid, 'purple', linewidth=0.5)
axes[2].set_title(f'LNK-Lipid Contacts (<5Å): {cont_lnk_lipid.mean():.1f}±{cont_lnk_lipid.std():.1f}')
axes[2].set_ylabel('Contacts')
axes[2].set_xlabel('Time (ns)')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/contacts.png", dpi=150)
plt.close()
log(f"  Saved: {OUTDIR}/contacts.png")

# =============================================================================
# 10. PCA (Protein CA)
# =============================================================================
log("\n[10] PCA on Protein CA...")

try:
    # Collect aligned coordinates
    pca_coords = []
    for ts in u.trajectory[::STEP]:
        unwrap_group(sel_protein_ca)
        align.alignto(sel_protein_ca, ref_protein_ca, select='name CA', weights='mass')
        pca_coords.append(sel_protein_ca.positions.flatten())
    pca_coords = np.array(pca_coords)
    
    pca = MDA_PCA(pca_coords, n_components=10)
    variance = pca.results.variance
    cumvar = pca.results.cumulated_variance
    
    df_pca = pd.DataFrame({
        'PC': range(1, 11),
        'Variance': variance,
        'Cumulative_Variance': cumvar,
    })
    df_pca.to_csv(f"{OUTDIR}/pca_variance.csv", index=False)
    log(f"  Saved: {OUTDIR}/pca_variance.csv")
    
    # Project trajectory onto PC1/PC2
    projections = pca_coords.dot(pca.results.pca.components_[:2].T)
    
    df_pca_proj = pd.DataFrame({
        'Time_ns': times_pos[:len(projections)],
        'PC1': projections[:, 0],
        'PC2': projections[:, 1],
    })
    df_pca_proj.to_csv(f"{OUTDIR}/pca_projection.csv", index=False)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].bar(df_pca['PC'], df_pca['Variance'], color='steelblue')
    axes[0].set_xlabel('Principal Component')
    axes[0].set_ylabel('Variance')
    axes[0].set_title('PCA Variance (Protein CA)')
    
    scatter = axes[1].scatter(df_pca_proj['PC1'], df_pca_proj['PC2'], c=df_pca_proj['Time_ns'], cmap='viridis', s=5)
    axes[1].set_xlabel('PC1')
    axes[1].set_ylabel('PC2')
    axes[1].set_title('PCA Projection (colored by time)')
    plt.colorbar(scatter, ax=axes[1], label='Time (ns)')
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/pca_overview.png", dpi=150)
    plt.close()
    log(f"  Saved: {OUTDIR}/pca_overview.png")
    log(f"  PC1 variance: {variance[0]:.2f}, PC2 variance: {variance[1]:.2f}")
except Exception as e:
    log(f"  PCA failed: {e}")

# =============================================================================
# 11. SUMMARY REPORT
# =============================================================================
log("\n[11] Generating summary report...")

report = f"""
================================================================================
GROMACS MEMBRANE MD ANALYSIS SUMMARY
================================================================================
Trajectory: {XTC}
Frames analyzed: {n_analyzed} (step={STEP})
Simulation time: {times[-1]:.1f} ns

--- SYSTEM COMPOSITION ---
ECD residues:     {len(sel_ecd.residues)}
GLP-1 residues:   {len(sel_glp1.residues)}
LNK atoms:        {len(sel_lnk.atoms)}
Lipid residues:   {len(sel_lipid.residues)}
Water residues:   {len(sel_water.residues)}

--- RMSD (Aligned on ECD CA, Å) ---
Protein CA:       {df_rmsd['RMSD_Protein_CA_A'].mean():.2f} ± {df_rmsd['RMSD_Protein_CA_A'].std():.2f}
ECD CA:           {df_rmsd['RMSD_ECD_CA_A'].mean():.2f} ± {df_rmsd['RMSD_ECD_CA_A'].std():.2f}
GLP-1 CA:         {df_rmsd['RMSD_GLP1_CA_A'].mean():.2f} ± {df_rmsd['RMSD_GLP1_CA_A'].std():.2f}
LNK Heavy:        {df_rmsd['RMSD_LNK_Heavy_A'].mean():.2f} ± {df_rmsd['RMSD_LNK_Heavy_A'].std():.2f}

--- RADIUS OF GYRATION (Å) ---
Protein:          {df_rg['Rg_Protein_A'].mean():.1f} ± {df_rg['Rg_Protein_A'].std():.1f}
ECD:              {df_rg['Rg_ECD_A'].mean():.1f} ± {df_rg['Rg_ECD_A'].std():.1f}
GLP-1:            {df_rg['Rg_GLP1_A'].mean():.1f} ± {df_rg['Rg_GLP1_A'].std():.1f}

--- RMSF (Å) ---
ECD CA max:       {df_rmsf_ecd['RMSF_A'].max():.2f}  (resindex {df_rmsf_ecd['RMSF_A'].idxmax()})
GLP-1 CA max:     {df_rmsf_glp1['RMSF_A'].max():.2f}  (resindex {df_rmsf_glp1['RMSF_A'].idxmax()})

--- LNK POSITION ---
Relative to mem:  {(lnk_com_z - mem_com_z).mean()*10:.2f} ± {(lnk_com_z - mem_com_z).std()*10:.2f} Å
Tilt angle:       {lnk_tilt.mean():.1f} ± {lnk_tilt.std():.1f}°
"""

if energy_ok and 'Box-X' in df_energy.columns:
    report += f"""
--- MEMBRANE PROPERTIES ---
Area per lipid:   {df_mem['Area_per_lipid_A2'].mean():.1f} ± {df_mem['Area_per_lipid_A2'].std():.1f} Å²
Thickness:        {membrane_thickness.mean():.1f} ± {membrane_thickness.std():.1f} Å
"""

if hb_internal is not None:
    report += f"""
--- H-BONDS ---
Protein internal: {df_hb['HB_Internal'].mean():.1f} ± {df_hb['HB_Internal'].std():.1f}
Protein-Lipid:    {df_hb['HB_Protein_Lipid'].mean():.1f} ± {df_hb['HB_Protein_Lipid'].std():.1f}
"""

report += f"""
--- CONTACTS (<5Å) ---
ECD-Lipid:        {cont_ecd_lipid.mean():.1f} ± {cont_ecd_lipid.std():.1f}
GLP-1-Lipid:      {cont_glp1_lipid.mean():.1f} ± {cont_glp1_lipid.std():.1f}
LNK-Lipid:        {cont_lnk_lipid.mean():.1f} ± {cont_lnk_lipid.std():.1f}
"""

if 'pca_variance.csv' in os.listdir(OUTDIR):
    report += f"""
--- PCA (Protein CA) ---
PC1 variance:     {variance[0]:.4f}
PC2 variance:     {variance[1]:.4f}
Top 5 cumulative: {cumvar[4]:.4f}
"""

report += """
================================================================================
"""

with open(f"{OUTDIR}/summary.txt", "w") as f:
    f.write(report)
log(report)

log("\n" + "="*70)
log("ANALYSIS COMPLETE")
log("="*70)
log(f"All outputs saved to: {OUTDIR}")
log_file.close()
