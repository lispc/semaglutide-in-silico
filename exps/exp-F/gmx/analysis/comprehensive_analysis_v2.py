#!/usr/bin/env python3
"""
Comprehensive GROMACS membrane MD analysis script v2.
Uses PBC-corrected protein trajectory for protein analysis.

Usage:
    cd exps/exp-F/gmx && python analysis/comprehensive_analysis_v2.py

Selections:
  ECD:      resindex 0-1112   (1113 residues)
  GLP-1:    resindex 1113-1138 (26 residues)
  LNK:      resindex 1139     (1 residue)
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.analysis.pca import PCA as MDA_PCA

warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)

# =============================================================================
# CONFIGURATION
# =============================================================================
REPO = "/home/scroll/personal/semaglutide-in-silico"
GMX_DIR = f"{REPO}/exps/exp-F/gmx"
OUTDIR = f"{GMX_DIR}/analysis"

# Full system (for membrane analysis)
GRO_FULL = f"{GMX_DIR}/npt.gro"
XTC_FULL = f"{GMX_DIR}/md.xtc"
EDR = f"{GMX_DIR}/md.edr"

# PBC-corrected protein trajectory (for protein analysis)
GRO_PROT = f"{OUTDIR}/protein_ref.gro"
XTC_PROT = f"{OUTDIR}/md_protein_whole.xtc"

STEP = 10  # every 10 frames = 20 ps

os.makedirs(OUTDIR, exist_ok=True)

log_file = open(f"{OUTDIR}/analysis_v2.log", "w")

def log(msg):
    print(msg)
    log_file.write(msg + "\n")
    log_file.flush()

log("="*70)
log("COMPREHENSIVE GROMACS MEMBRANE MD ANALYSIS v2")
log("="*70)

# =============================================================================
# 1. LOAD UNIVERSES
# =============================================================================
log("\n[1] Loading universes...")

u_full = mda.Universe(GRO_FULL, XTC_FULL)
n_frames_full = len(u_full.trajectory)
log(f"  Full system:   {len(u_full.atoms):6d} atoms, {n_frames_full:5d} frames")

u_prot = mda.Universe(GRO_PROT, XTC_PROT)
n_frames_prot = len(u_prot.trajectory)
log(f"  Protein whole: {len(u_prot.atoms):6d} atoms, {n_frames_prot:5d} frames")

# Reference structures
ref_full = mda.Universe(GRO_FULL)
ref_prot = mda.Universe(GRO_PROT)

# --- Full system selections ---
sel_lipid_full = u_full.select_atoms('resname PA PC OL CHL')
sel_lnk_full   = u_full.select_atoms('resindex 1139')
sel_ecd_full   = u_full.select_atoms('resindex 0-1112')
sel_glp1_full  = u_full.select_atoms('resindex 1113-1138')

# --- Protein selections (PBC-corrected) ---
sel_prot     = u_prot.select_atoms('protein')
sel_ecd      = u_prot.select_atoms('resindex 0-1112')
sel_glp1     = u_prot.select_atoms('resindex 1113-1138')
sel_lnk      = u_prot.select_atoms('resindex 1139')

sel_prot_ca  = sel_prot.select_atoms('name CA')
sel_ecd_ca   = sel_ecd.select_atoms('name CA')
sel_glp1_ca  = sel_glp1.select_atoms('name CA')

ref_prot_ca  = ref_prot.select_atoms('protein and name CA')
ref_ecd_ca   = ref_prot.select_atoms('resindex 0-1112 and name CA')
ref_glp1_ca  = ref_prot.select_atoms('resindex 1113-1138 and name CA')
ref_lnk      = ref_prot.select_atoms('resindex 1139')

log(f"  ECD:    {len(sel_ecd.atoms):6d} atoms, {len(sel_ecd.residues):4d} res")
log(f"  GLP-1:  {len(sel_glp1.atoms):6d} atoms, {len(sel_glp1.residues):4d} res")
log(f"  LNK:    {len(sel_lnk.atoms):6d} atoms, {len(sel_lnk.residues):4d} res")
log(f"  Lipid:  {len(sel_lipid_full.atoms):6d} atoms, {len(sel_lipid_full.residues):4d} res")

# =============================================================================
# 2. ENERGY
# =============================================================================
log("\n[2] Energy analysis...")
try:
    import panedr
    df_energy = panedr.edr_to_df(EDR)
    avail = [c for c in ['Potential','Temperature','Pressure','Volume','Density','Box-X','Box-Y','Box-Z'] if c in df_energy.columns]
    df_energy[['Time'] + avail].to_csv(f"{OUTDIR}/energy.csv", index=False)
    log(f"  Saved: {OUTDIR}/energy.csv")
    
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    t_ns = df_energy['Time'].values / 1000.0
    for ax, (col, label) in zip(axes.flat, [
        ('Potential', 'Potential (kJ/mol)'),
        ('Temperature', 'Temperature (K)'),
        ('Pressure', 'Pressure (bar)'),
        ('Density', 'Density (kg/m³)')
    ]):
        if col in df_energy.columns:
            ax.plot(t_ns, df_energy[col].values, linewidth=0.5)
            ax.set_title(label)
            ax.set_xlabel('Time (ns)')
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/energy_overview.png", dpi=150)
    plt.close()
    log(f"  Saved: energy_overview.png")
    energy_ok = True
except Exception as e:
    log(f"  Energy failed: {e}")
    energy_ok = False

# =============================================================================
# 3. RMSD
# =============================================================================
log("\n[3] RMSD analysis (aligned on ECD CA)...")

def aligned_rmsd(u_obj, ref_obj, ag, ref_ag, label, align_sel='name CA', step=STEP):
    vals = []
    times = []
    ref_pos = ref_ag.positions.copy()
    for ts in u_obj.trajectory[::step]:
        if align_sel:
            align.alignto(ag, ref_ag, select=align_sel, weights='mass')
        else:
            align.alignto(ag, ref_ag)
        vals.append(rms.rmsd(ag.positions, ref_pos, superposition=False))
        times.append(ts.time / 1000.0)
    return np.array(times), np.array(vals)

t_rmsd, rmsd_prot = aligned_rmsd(u_prot, ref_prot, sel_prot_ca, ref_prot_ca, 'Protein')
t_rmsd, rmsd_ecd  = aligned_rmsd(u_prot, ref_prot, sel_ecd_ca,  ref_ecd_ca,  'ECD')
t_rmsd, rmsd_glp1 = aligned_rmsd(u_prot, ref_prot, sel_glp1_ca, ref_glp1_ca, 'GLP-1')
t_rmsd, rmsd_lnk  = aligned_rmsd(u_prot, ref_prot, sel_lnk,     ref_lnk,     'LNK', align_sel=None)

df_rmsd = pd.DataFrame({
    'Time_ns': t_rmsd,
    'Protein_CA_A': rmsd_prot,
    'ECD_CA_A': rmsd_ecd,
    'GLP1_CA_A': rmsd_glp1,
    'LNK_Heavy_A': rmsd_lnk,
})
df_rmsd.to_csv(f"{OUTDIR}/rmsd.csv", index=False)

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
for ax, (col, title) in zip(axes.flat, [
    ('Protein_CA_A', 'Protein CA'),
    ('ECD_CA_A', 'ECD CA'),
    ('GLP1_CA_A', 'GLP-1 CA'),
    ('LNK_Heavy_A', 'LNK Heavy')
]):
    ax.plot(df_rmsd['Time_ns'], df_rmsd[col], linewidth=0.5)
    ax.set_title(f'{title}: {df_rmsd[col].mean():.1f}±{df_rmsd[col].std():.1f} Å')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsd_overview.png", dpi=150)
plt.close()
log(f"  Saved: rmsd.csv, rmsd_overview.png")

# =============================================================================
# 4. RADIUS OF GYRATION
# =============================================================================
log("\n[4] Radius of Gyration...")

def compute_rg(u_obj, selection, step=STEP):
    rg_vals, times = [], []
    for ts in u_obj.trajectory[::step]:
        rg_vals.append(selection.radius_of_gyration())
        times.append(ts.time / 1000.0)
    return np.array(times), np.array(rg_vals)

t_rg, rg_prot = compute_rg(u_prot, sel_prot, STEP)
_, rg_ecd = compute_rg(u_prot, sel_ecd, STEP)
_, rg_glp1 = compute_rg(u_prot, sel_glp1, STEP)

df_rg = pd.DataFrame({
    'Time_ns': t_rg,
    'Protein_A': rg_prot,
    'ECD_A': rg_ecd,
    'GLP1_A': rg_glp1,
})
df_rg.to_csv(f"{OUTDIR}/rg.csv", index=False)

fig, axes = plt.subplots(1, 3, figsize=(12, 4))
for ax, (col, title) in zip(axes, [
    ('Protein_A', 'Protein'), ('ECD_A', 'ECD'), ('GLP1_A', 'GLP-1')
]):
    ax.plot(df_rg['Time_ns'], df_rg[col], linewidth=0.5)
    ax.set_title(f'{title}: {df_rg[col].mean():.1f}±{df_rg[col].std():.1f} Å')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Rg (Å)')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rg_overview.png", dpi=150)
plt.close()
log(f"  Saved: rg.csv, rg_overview.png")

# =============================================================================
# 5. RMSF
# =============================================================================
log("\n[5] RMSF analysis...")

def compute_rmsf(u_obj, ag, ref_ag, step=STEP):
    coords = []
    for ts in u_obj.trajectory[::step]:
        align.alignto(ag, ref_ag, select='name CA', weights='mass')
        coords.append(ag.positions.copy())
    coords = np.array(coords)
    mean_pos = coords.mean(axis=0)
    return np.sqrt(((coords - mean_pos)**2).sum(axis=2)).mean(axis=0)

rmsf_ecd = compute_rmsf(u_prot, sel_ecd_ca, ref_ecd_ca, STEP)
rmsf_glp1 = compute_rmsf(u_prot, sel_glp1_ca, ref_glp1_ca, STEP)

df_rmsf_ecd = pd.DataFrame({
    'ResIndex': range(len(sel_ecd_ca)),
    'ResID': [r.resid for r in sel_ecd_ca.residues],
    'ResName': [r.resname for r in sel_ecd_ca.residues],
    'RMSF_A': rmsf_ecd,
})
df_rmsf_ecd.to_csv(f"{OUTDIR}/rmsf_ecd.csv", index=False)

df_rmsf_glp1 = pd.DataFrame({
    'ResIndex': range(len(sel_glp1_ca)),
    'ResID': [r.resid for r in sel_glp1_ca.residues],
    'ResName': [r.resname for r in sel_glp1_ca.residues],
    'RMSF_A': rmsf_glp1,
})
df_rmsf_glp1.to_csv(f"{OUTDIR}/rmsf_glp1.csv", index=False)

fig, axes = plt.subplots(2, 1, figsize=(12, 8))
axes[0].plot(df_rmsf_ecd['ResIndex'], df_rmsf_ecd['RMSF_A'], 'b-', linewidth=0.5)
axes[0].set_title('ECD CA RMSF (Å)')
axes[0].set_ylabel('RMSF (Å)')

axes[1].plot(df_rmsf_glp1['ResIndex'], df_rmsf_glp1['RMSF_A'], 'g-o', linewidth=1, markersize=3)
for i, (ri, rn, val) in enumerate(zip(df_rmsf_glp1['ResIndex'], df_rmsf_glp1['ResName'], df_rmsf_glp1['RMSF_A'])):
    if val > df_rmsf_glp1['RMSF_A'].mean() + 1.5 * df_rmsf_glp1['RMSF_A'].std():
        axes[1].annotate(f'{rn}{ri+1}', xy=(ri, val), xytext=(ri, val+0.5), fontsize=7, ha='center')
axes[1].set_title('GLP-1 CA RMSF (Å)')
axes[1].set_ylabel('RMSF (Å)')
axes[1].set_xlabel('Residue')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/rmsf_overview.png", dpi=150)
plt.close()
log(f"  Saved: rmsf_ecd.csv, rmsf_glp1.csv, rmsf_overview.png")

# =============================================================================
# 6. LNK POSITION (using full trajectory)
# =============================================================================
log("\n[6] LNK position & membrane analysis...")

phosphorus = u_full.select_atoms('name P31 or name P')
if len(phosphorus) == 0:
    phosphorus = u_full.select_atoms('name P* and (resname PA PC OL)')

lnk_head = sel_lnk_full.select_atoms('name C13 C14 C15 N16')
lnk_tail = sel_lnk_full.select_atoms('name C48 C49 C50 C51 C52 C53 C54 C55')
lys_nz   = u_full.select_atoms('resindex 1129 and name NZ')

lnk_com_z, lnk_head_z, lnk_tail_z = [], [], []
mem_upper_z, mem_lower_z, mem_com_z = [], [], []
lnk_tilt, lys_lnk_dist = [], []
times_pos = []

for ts in u_full.trajectory[::STEP]:
    p_z = phosphorus.positions[:, 2]
    upper, lower = np.percentile(p_z, 90), np.percentile(p_z, 10)
    mem_com = p_z.mean()
    
    lnk_com = sel_lnk_full.center_of_mass()
    lh = lnk_head.center_of_mass() if len(lnk_head) > 0 else lnk_com
    lt = lnk_tail.center_of_mass() if len(lnk_tail) > 0 else lnk_com
    
    vec = lt - lh
    vec_n = np.linalg.norm(vec)
    tilt = np.degrees(np.arccos(abs(vec[2]) / vec_n)) if vec_n > 0 else 0
    
    d = np.linalg.norm(lys_nz.center_of_mass() - lnk_com) if len(lys_nz) > 0 else np.nan
    
    lnk_com_z.append(lnk_com[2])
    lnk_head_z.append(lh[2])
    lnk_tail_z.append(lt[2])
    mem_upper_z.append(upper)
    mem_lower_z.append(lower)
    mem_com_z.append(mem_com)
    lnk_tilt.append(tilt)
    lys_lnk_dist.append(d)
    times_pos.append(ts.time / 1000.0)

times_pos = np.array(times_pos)
df_pos = pd.DataFrame({
    'Time_ns': times_pos,
    'LNK_COM_z_A': np.array(lnk_com_z),
    'LNK_Head_z_A': np.array(lnk_head_z),
    'LNK_Tail_z_A': np.array(lnk_tail_z),
    'Mem_Upper_z_A': np.array(mem_upper_z),
    'Mem_Lower_z_A': np.array(mem_lower_z),
    'Mem_COM_z_A': np.array(mem_com_z),
    'LNK_Tilt_deg': lnk_tilt,
    'LysNZ_LNK_Dist_A': np.array(lys_lnk_dist),
    'LNK_COM_rel_mem_A': (np.array(lnk_com_z) - np.array(mem_com_z)),
})
df_pos.to_csv(f"{OUTDIR}/lnk_position.csv", index=False)

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
axes[0,0].plot(times_pos, df_pos['LNK_COM_z_A'], 'purple', linewidth=0.5, label='LNK COM')
axes[0,0].plot(times_pos, df_pos['Mem_Upper_z_A'], 'b--', alpha=0.5, linewidth=0.5, label='Mem upper')
axes[0,0].plot(times_pos, df_pos['Mem_Lower_z_A'], 'b--', alpha=0.5, linewidth=0.5, label='Mem lower')
axes[0,0].set_title('LNK Z Position (Å)')
axes[0,0].legend(fontsize=7)

axes[0,1].plot(times_pos, df_pos['LNK_COM_rel_mem_A'], 'purple', linewidth=0.5)
axes[0,1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
axes[0,1].set_title('LNK COM relative to Membrane COM (Å)')

axes[1,0].plot(times_pos, df_pos['LNK_Tilt_deg'], 'red', linewidth=0.5)
axes[1,0].set_title('LNK Tilt Angle (°)')
axes[1,0].set_xlabel('Time (ns)')

axes[1,1].plot(times_pos, df_pos['LysNZ_LNK_Dist_A'], 'green', linewidth=0.5)
axes[1,1].set_title('Lys NZ — LNK Distance (Å)')
axes[1,1].set_xlabel('Time (ns)')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/lnk_position.png", dpi=150)
plt.close()
log(f"  Saved: lnk_position.csv, lnk_position.png")

# =============================================================================
# 7. MEMBRANE PROPERTIES
# =============================================================================
log("\n[7] Membrane properties...")

mem_thickness = (np.array(mem_upper_z) - np.array(mem_lower_z))
n_lipids = len(phosphorus.residues) if len(phosphorus) > 0 else len(sel_lipid_full.residues)

if energy_ok and 'Box-X' in df_energy.columns:
    box_x = df_energy['Box-X'].values
    box_y = df_energy['Box-Y'].values if 'Box-Y' in df_energy.columns else box_x
    area_per_lipid = (box_x * box_y * 100) / (n_lipids / 2)
    
    df_mem = pd.DataFrame({
        'Time_ns': df_energy['Time'].values / 1000.0,
        'Area_per_lipid_A2': area_per_lipid,
        'Membrane_thickness_A': np.interp(df_energy['Time'].values / 1000.0, times_pos, mem_thickness),
    })
    df_mem.to_csv(f"{OUTDIR}/membrane_properties.csv", index=False)
    
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    axes[0].plot(df_mem['Time_ns'], df_mem['Area_per_lipid_A2'], 'b-', linewidth=0.5)
    axes[0].set_title(f'Area/Lipid: {df_mem["Area_per_lipid_A2"].mean():.1f}±{df_mem["Area_per_lipid_A2"].std():.1f} Å²')
    axes[1].plot(times_pos, mem_thickness, 'b-', linewidth=0.5)
    axes[1].set_title(f'Thickness: {mem_thickness.mean():.1f}±{mem_thickness.std():.1f} Å')
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/membrane_properties.png", dpi=150)
    plt.close()
    log(f"  Saved: membrane_properties.csv, membrane_properties.png")
else:
    log("  Skipped (no box data)")

# =============================================================================
# 8. CONTACTS (fast: CA vs lipid P atoms)
# =============================================================================
log("\n[8] Protein-lipid contacts (CA vs P headgroups)...")

# Select lipid phosphorus atoms for headgroup contact analysis
lipid_p = u_full.select_atoms('name P31 or name P')
if len(lipid_p) == 0:
    lipid_p = u_full.select_atoms('name P* and (resname PA PC OL)')

# Protein CA selections for contact analysis
ecd_ca_full   = sel_ecd_full.select_atoms('name CA')
glp1_ca_full  = sel_glp1_full.select_atoms('name CA')
lnk_full      = sel_lnk_full  # LNK has no CA, use all heavy atoms

def compute_contacts(u_obj, sel1, sel2, cutoff_A=5.0, step=STEP):
    log(f"  Contacts: {len(sel1.atoms)} atoms - {len(sel2.atoms)} atoms...")
    counts, times_c = [], []
    for ts in u_obj.trajectory[::step]:
        pairs = mda.lib.distances.capped_distance(
            sel1.positions, sel2.positions,
            max_cutoff=cutoff_A,
            return_distances=False,
            box=u_obj.dimensions,
        )
        counts.append(len(pairs))
        times_c.append(ts.time / 1000.0)
    return np.array(times_c), np.array(counts)

t_c, cont_ecd = compute_contacts(u_full, ecd_ca_full, lipid_p, 5.0, STEP)
_, cont_glp1 = compute_contacts(u_full, glp1_ca_full, lipid_p, 5.0, STEP)
_, cont_lnk = compute_contacts(u_full, lnk_full, lipid_p, 5.0, STEP)

df_cont = pd.DataFrame({
    'Time_ns': t_c,
    'Contacts_ECD_Lipid': cont_ecd,
    'Contacts_GLP1_Lipid': cont_glp1,
    'Contacts_LNK_Lipid': cont_lnk,
})
df_cont.to_csv(f"{OUTDIR}/contacts.csv", index=False)

fig, axes = plt.subplots(3, 1, figsize=(10, 8))
for ax, (col, title) in zip(axes, [
    ('Contacts_ECD_Lipid', 'ECD CA - Lipid P'),
    ('Contacts_GLP1_Lipid', 'GLP-1 CA - Lipid P'),
    ('Contacts_LNK_Lipid', 'LNK Heavy - Lipid P')
]):
    ax.plot(t_c, df_cont[col], linewidth=0.5)
    ax.set_title(f'{title}: {df_cont[col].mean():.1f}±{df_cont[col].std():.1f}')
    ax.set_ylabel('Contacts (<5Å)')
axes[-1].set_xlabel('Time (ns)')
plt.tight_layout()
plt.savefig(f"{OUTDIR}/contacts.png", dpi=150)
plt.close()
log(f"  Saved: contacts.csv, contacts.png")

# =============================================================================
# 9. H-BONDS (lightweight: only GLP-1-ECD interface)
# =============================================================================
log("\n[9] Hydrogen bond analysis (GLP-1-ECD interface only)...")

def run_hbonds(u_obj, sel1_str, sel2_str, label, step=STEP):
    log(f"  {label}...")
    try:
        hb = HydrogenBondAnalysis(
            universe=u_obj,
            donors_sel=sel1_str,
            hydrogens_sel=sel1_str,
            acceptors_sel=sel2_str,
            d_a_cutoff=3.5,
            d_h_a_angle_cutoff=150,
            update_selections=False,
        )
        hb.run(step=step, verbose=False)
        n_frames = len(u_obj.trajectory[::step])
        counts = [0] * n_frames
        if len(hb.results.hbonds) > 0:
            for row in hb.results.hbonds:
                counts[int(row[0])] += 1
        return np.array(counts)
    except Exception as e:
        log(f"    Failed: {e}")
        return None

try:
    hb_glp1_ecd = run_hbonds(u_prot, 'resindex 1113-1138', 'resindex 0-1112', 'GLP-1-ECD', STEP)
except Exception as e:
    log(f"  H-bond analysis failed: {e}")
    hb_glp1_ecd = None

if hb_glp1_ecd is not None:
    df_hb = pd.DataFrame({
        'Time_ns': t_rmsd[:len(hb_glp1_ecd)],
        'HB_GLP1_ECD': hb_glp1_ecd,
    })
    df_hb.to_csv(f"{OUTDIR}/hbonds.csv", index=False)
    
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(df_hb['Time_ns'], df_hb['HB_GLP1_ECD'], linewidth=0.5)
    ax.set_title(f'GLP-1-ECD H-Bonds: {df_hb["HB_GLP1_ECD"].mean():.1f}±{df_hb["HB_GLP1_ECD"].std():.1f}')
    ax.set_ylabel('Count')
    ax.set_xlabel('Time (ns)')
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/hbonds.png", dpi=150)
    plt.close()
    log(f"  Saved: hbonds.csv, hbonds.png")

# =============================================================================
# 10. PCA
# =============================================================================
log("\n[10] PCA on Protein CA...")
try:
    pca_coords = []
    for ts in u_prot.trajectory[::STEP]:
        align.alignto(sel_prot_ca, ref_prot_ca, select='name CA', weights='mass')
        pca_coords.append(sel_prot_ca.positions.flatten())
    pca_coords = np.array(pca_coords)
    
    # Manual PCA via SVD (no sklearn dependency)
    pca_coords_centered = pca_coords - pca_coords.mean(axis=0)
    U, s, Vt = np.linalg.svd(pca_coords_centered, full_matrices=False)
    var = (s**2) / np.sum(s**2)
    cumvar = np.cumsum(var)
    
    df_pca_var = pd.DataFrame({'PC': range(1, 11), 'Variance': var, 'Cumulative': cumvar})
    df_pca_var.to_csv(f"{OUTDIR}/pca_variance.csv", index=False)
    
    # Project onto first 2 PCs
    proj = pca_coords_centered.dot(Vt[:2].T)
    df_pca_proj = pd.DataFrame({'Time_ns': t_rmsd[:len(proj)], 'PC1': proj[:, 0], 'PC2': proj[:, 1]})
    df_pca_proj.to_csv(f"{OUTDIR}/pca_projection.csv", index=False)
    
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    axes[0].bar(df_pca_var['PC'], df_pca_var['Variance'], color='steelblue')
    axes[0].set_xlabel('PC')
    axes[0].set_ylabel('Variance')
    axes[0].set_title('PCA Variance')
    
    sc = axes[1].scatter(df_pca_proj['PC1'], df_pca_proj['PC2'], c=df_pca_proj['Time_ns'], cmap='viridis', s=5)
    axes[1].set_xlabel('PC1')
    axes[1].set_ylabel('PC2')
    axes[1].set_title('PCA Projection')
    plt.colorbar(sc, ax=axes[1], label='Time (ns)')
    plt.tight_layout()
    plt.savefig(f"{OUTDIR}/pca_overview.png", dpi=150)
    plt.close()
    log(f"  Saved: pca_variance.csv, pca_projection.csv, pca_overview.png")
    log(f"  PC1: {var[0]:.4f}, PC2: {var[1]:.4f}, Top5 cum: {cumvar[4]:.4f}")
except Exception as e:
    log(f"  PCA failed: {e}")

# =============================================================================
# 11. SUMMARY
# =============================================================================
log("\n[11] Summary report...")

report = f"""================================================================================
GROMACS MEMBRANE MD ANALYSIS SUMMARY
================================================================================
Trajectory: {XTC_FULL}
Protein trajectory: {XTC_PROT}
Frames analyzed: {len(t_rmsd)} (step={STEP})
Simulation time: {t_rmsd[-1]:.1f} ns

--- RMSD (Aligned on ECD CA, Å) ---
Protein CA:  {df_rmsd['Protein_CA_A'].mean():.2f} ± {df_rmsd['Protein_CA_A'].std():.2f}
ECD CA:      {df_rmsd['ECD_CA_A'].mean():.2f} ± {df_rmsd['ECD_CA_A'].std():.2f}
GLP-1 CA:    {df_rmsd['GLP1_CA_A'].mean():.2f} ± {df_rmsd['GLP1_CA_A'].std():.2f}
LNK Heavy:   {df_rmsd['LNK_Heavy_A'].mean():.2f} ± {df_rmsd['LNK_Heavy_A'].std():.2f}

--- RADIUS OF GYRATION (Å) ---
Protein:     {df_rg['Protein_A'].mean():.1f} ± {df_rg['Protein_A'].std():.1f}
ECD:         {df_rg['ECD_A'].mean():.1f} ± {df_rg['ECD_A'].std():.1f}
GLP-1:       {df_rg['GLP1_A'].mean():.1f} ± {df_rg['GLP1_A'].std():.1f}

--- RMSF (Å) ---
ECD CA max:  {df_rmsf_ecd['RMSF_A'].max():.2f} (resindex {df_rmsf_ecd['RMSF_A'].idxmax()})
GLP-1 CA max: {df_rmsf_glp1['RMSF_A'].max():.2f} (resindex {df_rmsf_glp1['RMSF_A'].idxmax()})

--- LNK POSITION ---
Relative to mem: {df_pos['LNK_COM_rel_mem_A'].mean():.2f} ± {df_pos['LNK_COM_rel_mem_A'].std():.2f} Å
Tilt angle:      {df_pos['LNK_Tilt_deg'].mean():.1f} ± {df_pos['LNK_Tilt_deg'].std():.1f}°
"""

if energy_ok and 'Box-X' in df_energy.columns:
    report += f"""
--- MEMBRANE ---
Area/lipid:  {df_mem['Area_per_lipid_A2'].mean():.1f} ± {df_mem['Area_per_lipid_A2'].std():.1f} Å²
Thickness:   {mem_thickness.mean():.1f} ± {mem_thickness.std():.1f} Å
"""

if 'df_hb' in locals() and hb_glp1_ecd is not None:
    report += f"""
--- H-BONDS ---
GLP-1-ECD:   {df_hb['HB_GLP1_ECD'].mean():.1f} ± {df_hb['HB_GLP1_ECD'].std():.1f}
"""

report += f"""
--- CONTACTS (<5Å) ---
ECD-Lipid:   {cont_ecd.mean():.1f} ± {cont_ecd.std():.1f}
GLP-1-Lipid: {cont_glp1.mean():.1f} ± {cont_glp1.std():.1f}
LNK-Lipid:   {cont_lnk.mean():.1f} ± {cont_lnk.std():.1f}
"""

if os.path.exists(f"{OUTDIR}/pca_variance.csv"):
    report += f"""
--- PCA ---
PC1 variance: {var[0]:.4f}
PC2 variance: {var[1]:.4f}
Top5 cum:     {cumvar[4]:.4f}
"""

report += "="*70 + "\n"

with open(f"{OUTDIR}/summary.txt", "w") as f:
    f.write(report)
log(report)

log("\n" + "="*70)
log("ANALYSIS COMPLETE")
log("="*70)
log(f"All outputs: {OUTDIR}")
log_file.close()
