#!/usr/bin/env python3
"""
Analyze ECD v2 production MD trajectory.
Computes RMSD, RMSF, and thermodynamic metrics.

Usage:
    cd exps/exp-F/md/ecd_v2 && python ../../scripts/analyze_ecd_v2.py
"""
import sys, os
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.rms import RMSF
import numpy as np

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/minimal_model/system_ecd_v2.prmtop"
DCD = f"{EXP_F}/md/ecd_v2/prod.dcd"

print(f"Loading {DCD} ...")
u = mda.Universe(PRMTOP, DCD)
ref = u.copy()
ref.trajectory[0]
print(f"Frames: {len(u.trajectory)}")

# Components
hsa_ca = "name CA and resid 1-131"
pep_ca = "name CA and resid 132-145"
lnk_ca = "name CA and resid 146-157"

# ---------- Internal RMSD (self-aligned) ----------
def internal_rmsd(label, sel_str):
    u2 = mda.Universe(PRMTOP, DCD)
    ref2 = u2.copy()
    ref2.trajectory[0]
    align.AlignTraj(u2, ref2, select=sel_str, in_memory=True).run()
    r = rms.RMSD(u2, ref2, select=sel_str)
    r.run()
    vals = r.results.rmsd[:, 2]
    print(f"  {label:12s}: mean={np.mean(vals):.3f} Å  max={np.max(vals):.3f} Å  drift={np.mean(vals[-10:])-np.mean(vals[:10]):+.3f} Å")
    return vals

print("\n=== Internal RMSD (self-aligned) ===")
internal_rmsd("HSA", hsa_ca)
internal_rmsd("Peptide", pep_ca)
internal_rmsd("Linker", lnk_ca)

# ---------- Internal RMSF (self-aligned) ----------
def internal_rmsf(label, sel_str):
    u2 = mda.Universe(PRMTOP, DCD)
    ref2 = u2.copy()
    ref2.trajectory[0]
    align.AlignTraj(u2, ref2, select=sel_str, in_memory=True).run()
    ca = u2.select_atoms(sel_str)
    r = RMSF(ca).run()
    vals = r.results.rmsf
    print(f"  {label:12s}: mean={np.mean(vals):.3f} Å  max={np.max(vals):.3f} Å")
    return vals, ca.resids, ca.resnames

print("\n=== Internal RMSF (self-aligned) ===")
internal_rmsf("HSA", hsa_ca)
internal_rmsf("Peptide", pep_ca)
internal_rmsf("Linker", lnk_ca)

# ---------- Thermodynamic log analysis ----------
print("\n=== Production Log Metrics ===")
log_path = f"{EXP_F}/md/ecd_v2/prod.log"
if os.path.exists(log_path):
    import pandas as pd
    prod = pd.read_csv(log_path, skiprows=1, header=None)
    prod.columns = ['Step','Time_ps','PE','KE','TE','Temp','Vol','Dens','Speed']
    print(f"  Frames: {len(prod)}")
    print(f"  Time: {prod['Time_ps'].iloc[0]/1000:.2f} - {prod['Time_ps'].iloc[-1]/1000:.2f} ns")
    print(f"  Temp: {prod['Temp'].mean():.2f} ± {prod['Temp'].std():.2f} K")
    print(f"  Density: {prod['Dens'].mean():.4f} ± {prod['Dens'].std():.4f} g/mL")
    print(f"  Speed: {prod['Speed'].iloc[-1]:.1f} ns/day")
    remaining_ns = (50_000_000 - prod['Step'].iloc[-1]) / 500_000
    remaining_h = remaining_ns / prod['Speed'].iloc[-1] * 24
    print(f"  Progress: {prod['Step'].iloc[-1]/50_000_000*100:.1f}%")
    print(f"  Remaining: ~{remaining_h:.1f} h")
