#!/usr/bin/env python3
"""Analyze solvent and membrane MD trajectories for structural stability."""
import sys, os
import numpy as np
import mdtraj as md

REPO = "/home/scroll/personal/semaglutide-in-silico"

def analyze_system(name, prmtop_path, dcd_path, is_membrane=False):
    print(f"\n{'='*60}")
    print(f"ANALYZING: {name}")
    print(f"{'='*60}")
    
    # Load topology
    top = md.load_prmtop(prmtop_path)
    
    # Find key atom indices in original topology
    # Receptor CA: residues 0-1112
    receptor_ca = [a.index for a in top.atoms if a.name == 'CA' and a.residue.index < 1113]
    # Peptide CA: residues 1113-1138 (excluding LNK)
    peptide_ca = [a.index for a in top.atoms if a.name == 'CA' and 1113 <= a.residue.index <= 1138]
    # LYS NZ (residue 1129)
    lys_nz = [a.index for a in top.atoms if a.name == 'NZ' and a.residue.index == 1129]
    # LNK C13 (residue 1139)
    lnk_c13 = [a.index for a in top.atoms if a.name == 'C13' and a.residue.index == 1139]
    # Peptide N-terminus N (residue 1113)
    pep_n = [a.index for a in top.atoms if a.name == 'N' and a.residue.index == 1113]
    
    print(f"Receptor CA atoms: {len(receptor_ca)}")
    print(f"Peptide CA atoms: {len(peptide_ca)}")
    print(f"LYS NZ: {lys_nz}")
    print(f"LNK C13: {lnk_c13}")
    print(f"Peptide N: {pep_n}")
    
    # Create mapping from original index to position in loaded trajectory
    # We'll load all needed atoms at once
    all_needed = list(set(receptor_ca + peptide_ca + lys_nz + lnk_c13 + pep_n))
    
    # Load trajectory
    print(f"\nLoading trajectory ({len(all_needed)} atoms)...")
    traj_full = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=all_needed)
    print(f"Loaded {traj_full.n_frames} frames")
    
    # Build mapping: original atom index -> index in loaded trajectory
    orig_to_loaded = {orig_idx: loaded_idx for loaded_idx, orig_idx in enumerate(all_needed)}
    
    ref = traj_full[0]
    
    # Helper to get loaded indices
    def loaded_idx(orig_indices):
        return [orig_to_loaded[i] for i in orig_indices]
    
    # 1. Receptor CA RMSD
    rc_idx = loaded_idx(receptor_ca)
    rmsd_receptor = md.rmsd(traj_full.atom_slice(rc_idx), ref.atom_slice(rc_idx))
    print(f"\n[1] Receptor CA RMSD:")
    print(f"    Mean: {np.mean(rmsd_receptor):.3f} nm = {np.mean(rmsd_receptor)*10:.2f} Å")
    print(f"    Std:  {np.std(rmsd_receptor):.3f} nm")
    print(f"    Max:  {np.max(rmsd_receptor):.3f} nm = {np.max(rmsd_receptor)*10:.2f} Å")
    if np.max(rmsd_receptor) > 0.5:
        print(f"    ⚠️  High RMSD - receptor may be unstable")
    else:
        print(f"    ✅ Receptor stable")
    
    # 2. Peptide CA RMSD
    pc_idx = loaded_idx(peptide_ca)
    rmsd_peptide = md.rmsd(traj_full.atom_slice(pc_idx), ref.atom_slice(pc_idx))
    print(f"\n[2] Peptide CA RMSD:")
    print(f"    Mean: {np.mean(rmsd_peptide):.3f} nm = {np.mean(rmsd_peptide)*10:.2f} Å")
    print(f"    Std:  {np.std(rmsd_peptide):.3f} nm")
    print(f"    Max:  {np.max(rmsd_peptide):.3f} nm = {np.max(rmsd_peptide)*10:.2f} Å")
    if np.max(rmsd_peptide) > 0.5:
        print(f"    ⚠️  High RMSD - peptide conformation changing significantly")
    else:
        print(f"    ✅ Peptide conformation stable")
    
    # 3. LYS-LNK amide bond distance
    if lys_nz and lnk_c13:
        dist_bond = md.compute_distances(traj_full, [[orig_to_loaded[lys_nz[0]], orig_to_loaded[lnk_c13[0]]]])
        print(f"\n[3] LYS(1129) NZ - LNK(1139) C13 distance:")
        print(f"    Mean: {np.mean(dist_bond):.3f} nm = {np.mean(dist_bond)*10:.2f} Å")
        print(f"    Std:  {np.std(dist_bond):.3f} nm = {np.std(dist_bond)*10:.2f} Å")
        print(f"    Max:  {np.max(dist_bond):.3f} nm = {np.max(dist_bond)*10:.2f} Å")
        
        n_broken = np.sum(dist_bond > 0.25)
        if n_broken > 0:
            print(f"    ⚠️  WARNING: {n_broken}/{traj_full.n_frames} frames with bond > 2.5 Å")
        else:
            print(f"    ✅ Amide bond intact (all frames < 2.5 Å)")
    
    # 4. Peptide N-terminus to receptor distance (check for dissociation)
    if pep_n and receptor_ca:
        pn_idx = orig_to_loaded[pep_n[0]]
        rc_idx_in_full = loaded_idx(receptor_ca)
        min_dists = []
        for frame in range(traj_full.n_frames):
            pn_pos = traj_full.xyz[frame, pn_idx, :]
            ca_pos = traj_full.xyz[frame, rc_idx_in_full, :]
            dists = np.linalg.norm(ca_pos - pn_pos, axis=1)
            min_dists.append(np.min(dists))
        min_dists = np.array(min_dists)
        print(f"\n[4] Peptide N-terminus to receptor (min distance to any receptor CA):")
        print(f"    Mean: {np.mean(min_dists):.3f} nm = {np.mean(min_dists)*10:.2f} Å")
        print(f"    Std:  {np.std(min_dists):.3f} nm")
        print(f"    Initial: {min_dists[0]*10:.2f} Å")
        print(f"    Final:   {min_dists[-1]*10:.2f} Å")
        
        drift = min_dists[-1] - min_dists[0]
        if abs(drift) > 0.5:  # 5 Å
            print(f"    ⚠️  WARNING: N-terminus drifted {drift*10:.1f} Å - possible dissociation")
        else:
            print(f"    ✅ N-terminus position stable")
    
    # 5. Radius of gyration
    rg_peptide = md.compute_rg(traj_full.atom_slice(pc_idx))
    print(f"\n[5] Peptide Radius of Gyration:")
    print(f"    Mean: {np.mean(rg_peptide):.3f} nm")
    print(f"    Std:  {np.std(rg_peptide):.3f} nm")
    
    # 6. Membrane-specific analysis
    if is_membrane:
        print(f"\n[6] MEMBRANE ANALYSIS:")
        
        # Reload with lipid P atoms
        p_atoms = [a.index for a in top.atoms if a.name == 'P' and a.residue.index >= 1140]
        if p_atoms:
            print(f"    Lipid P atoms: {len(p_atoms)}")
            traj_p = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=p_atoms)
            
            # Bilayer thickness: z-distance between upper and lower leaflet P atoms
            z0 = traj_p.xyz[0, :, 2]
            mid_z = np.median(z0)
            upper_idx = np.where(z0 > mid_z)[0]
            lower_idx = np.where(z0 < mid_z)[0]
            
            thickness = []
            for frame in range(traj_p.n_frames):
                z = traj_p.xyz[frame, :, 2]
                upper_z = np.mean(z[upper_idx])
                lower_z = np.mean(z[lower_idx])
                thickness.append(abs(upper_z - lower_z))
            thickness = np.array(thickness)
            
            print(f"    Bilayer thickness (P-P distance):")
            print(f"      Mean: {np.mean(thickness):.3f} nm = {np.mean(thickness)*10:.1f} Å")
            print(f"      Std:  {np.std(thickness):.3f} nm")
            print(f"      Expected for POPC: ~3.7-4.0 nm")
            
            if np.mean(thickness) < 3.0 or np.mean(thickness) > 5.0:
                print(f"    ⚠️  WARNING: Bilayer thickness abnormal!")
            else:
                print(f"    ✅ Bilayer thickness normal")
        
        # LNK tail position
        if lnk_c13:
            lnk_atoms = [a for a in top.atoms if a.residue.index == 1139]
            tail_atoms = [a.index for a in lnk_atoms if a.name in ['C55', 'C58', 'O56', 'O57', 'O59']]
            if tail_atoms:
                traj_tail = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=tail_atoms)
                tail_z = traj_tail.xyz[:, :, 2]
                print(f"\n    LNK tail (diacid end) z-position:")
                print(f"      Mean: {np.mean(tail_z):.3f} nm")
                print(f"      Std:  {np.std(tail_z):.3f} nm")
                print(f"      Range: {np.min(tail_z):.3f} to {np.max(tail_z):.3f} nm")
                
                # Check if tail is in membrane or water
                if np.mean(np.abs(tail_z)) < 2.0:
                    print(f"    ⚠️  Tail near membrane center - indicates insertion into bilayer")
                elif np.mean(tail_z) > 3.5 or np.mean(tail_z) < -3.5:
                    print(f"    ✅ Tail in bulk aqueous phase (away from membrane)")
                else:
                    print(f"    ℹ️  Tail at membrane-water interface")

# Analyze solvent system
analyze_system(
    "Solvent System (exp-F Phase 1)",
    f"{REPO}/exps/exp-F/tleap/system.prmtop",
    f"{REPO}/exps/exp-F/md/rep1/prod.dcd",
    is_membrane=False
)

# Analyze membrane system
analyze_system(
    "Membrane System (exp-F Phase 2)",
    f"{REPO}/exps/exp-F/membrane_build/system_final.prmtop",
    f"{REPO}/exps/exp-F/md/membrane_rep1/prod.dcd",
    is_membrane=True
)
