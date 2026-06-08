#!/usr/bin/env python3
"""
Full structural analysis of solvent and membrane MD trajectories.
Treats current data as complete (69.8 ns solvent, 62.2 ns membrane).

Outputs: analysis/*.png, analysis/*.dat
"""
import sys, os, time
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

REPO = "/home/scroll/personal/semaglutide-in-silico"
OUT = f"{REPO}/exps/exp-F/analysis"
os.makedirs(OUT, exist_ok=True)

plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 9
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 11
plt.rcParams['legend.fontsize'] = 8
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8

def analyze(name, prmtop_path, dcd_path, is_membrane=False, max_frames=None):
    print(f"\n{'='*60}")
    print(f"FULL ANALYSIS: {name}")
    print(f"{'='*60}")
    t0 = time.time()
    
    top = md.load_prmtop(prmtop_path)
    
    # Atom selections
    receptor_ca = [a.index for a in top.atoms if a.name == 'CA' and a.residue.index < 1113]
    peptide_ca = [a.index for a in top.atoms if a.name == 'CA' and 1113 <= a.residue.index <= 1138]
    lnk_atoms = [a.index for a in top.atoms if a.residue.index == 1139]
    lys_nz = [a.index for a in top.atoms if a.name == 'NZ' and a.residue.index == 1129]
    lnk_c13 = [a.index for a in top.atoms if a.name == 'C13' and a.residue.index == 1139]
    pep_n = [a.index for a in top.atoms if a.name == 'N' and a.residue.index == 1113]
    
    # Determine stride for large trajectories (target ~500 frames for analysis)
    total_frames = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=[0]).n_frames
    stride = max(1, total_frames // 500)
    print(f"Total frames: {total_frames}, using stride={stride} for analysis (~{total_frames//stride} frames)")
    
    # === 1. RMSD analysis (CA only) ===
    print("[1/8] Computing RMSD...")
    traj_ca = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=receptor_ca + peptide_ca, stride=stride)
    ref = traj_ca[0]
    
    rc_idx = list(range(len(receptor_ca)))
    pc_idx = list(range(len(receptor_ca), len(receptor_ca) + len(peptide_ca)))
    
    rmsd_receptor = md.rmsd(traj_ca.atom_slice(rc_idx), ref.atom_slice(rc_idx))
    rmsd_peptide = md.rmsd(traj_ca.atom_slice(pc_idx), ref.atom_slice(pc_idx))
    
    # DCD timestamps from OpenMM are unreliable; derive time from frame index
    # Reporter saves every 5000 steps = 10 ps at 2 fs timestep
    time_ns = np.arange(traj_ca.n_frames) * 10.0 * stride / 1000.0  # ps -> ns
    
    fig, axes = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
    axes[0].plot(time_ns, rmsd_receptor * 10, color='tab:blue', linewidth=0.8)
    axes[0].axhline(y=np.mean(rmsd_receptor)*10, color='tab:blue', linestyle='--', alpha=0.5, label=f'Mean: {np.mean(rmsd_receptor)*10:.1f} Å')
    axes[0].set_ylabel('Receptor CA RMSD (Å)')
    axes[0].legend(loc='upper right')
    axes[0].set_title(f'{name}: Backbone Stability')
    
    axes[1].plot(time_ns, rmsd_peptide * 10, color='tab:red', linewidth=0.8)
    axes[1].axhline(y=np.mean(rmsd_peptide)*10, color='tab:red', linestyle='--', alpha=0.5, label=f'Mean: {np.mean(rmsd_peptide)*10:.1f} Å')
    axes[1].set_ylabel('Peptide CA RMSD (Å)')
    axes[1].set_xlabel('Time (ns)')
    axes[1].legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(f"{OUT}/{name.lower().replace(' ', '_')}_rmsd.png")
    plt.close()
    
    np.savetxt(f"{OUT}/{name.lower().replace(' ', '_')}_rmsd_receptor.dat", np.column_stack([time_ns, rmsd_receptor*10]), header='Time(ns) RMSD(A)')
    np.savetxt(f"{OUT}/{name.lower().replace(' ', '_')}_rmsd_peptide.dat", np.column_stack([time_ns, rmsd_peptide*10]), header='Time(ns) RMSD(A)')
    
    # === 2. RMSF analysis ===
    print("[2/8] Computing RMSF...")
    traj_prot = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=receptor_ca + peptide_ca, stride=stride)
    ref_prot = traj_prot[0]
    traj_prot.superpose(ref_prot, atom_indices=rc_idx)
    
    rmsf = md.rmsf(traj_prot, ref_prot, 0)
    res_indices = [traj_prot.topology.atom(i).residue.index for i in range(traj_prot.n_atoms)]
    
    # Separate receptor and peptide RMSF
    rec_rmsf = rmsf[:len(receptor_ca)] * 10
    pep_rmsf = rmsf[len(receptor_ca):] * 10
    
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(range(len(rec_rmsf)), rec_rmsf, color='tab:blue', linewidth=0.8, label='Receptor')
    ax.plot(range(len(rec_rmsf), len(rec_rmsf)+len(pep_rmsf)), pep_rmsf, color='tab:red', linewidth=0.8, label='Peptide')
    ax.axvline(x=len(rec_rmsf)-0.5, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Residue Index')
    ax.set_ylabel('RMSF (Å)')
    ax.set_title(f'{name}: Per-Residue Flexibility (RMSF)')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{OUT}/{name.lower().replace(' ', '_')}_rmsf.png")
    plt.close()
    
    # === 3. Key distances ===
    print("[3/8] Computing key distances...")
    
    # Use same trajectory (traj_ca) for distance calc to ensure frame alignment
    # Receptor CA indices in traj_ca: rc_idx, Peptide CA: pc_idx
    # Need to find LYS NZ, LNK C13, Peptide N in traj_ca topology
    
    # Build mapping from original atom index to traj_ca index
    all_ca_atoms = receptor_ca + peptide_ca
    ca_loaded_map = {orig: i for i, orig in enumerate(all_ca_atoms)}
    
    # LYS-LNK bond distance — need to load these specific atoms
    dist_extra_atoms = list(set(lys_nz + lnk_c13 + pep_n))
    traj_dist = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=dist_extra_atoms, stride=stride)
    dist_loaded_map = {orig: i for i, orig in enumerate(dist_extra_atoms)}
    
    if lys_nz and lnk_c13:
        dist_bond = md.compute_distances(traj_dist, [[dist_loaded_map[lys_nz[0]], dist_loaded_map[lnk_c13[0]]]])
        # Ensure same length as time_ns
        if len(dist_bond) > len(time_ns):
            dist_bond = dist_bond[:len(time_ns)]
        elif len(dist_bond) < len(time_ns):
            time_ns_dist = time_ns[:len(dist_bond)]
        else:
            time_ns_dist = time_ns
    else:
        dist_bond = None
        time_ns_dist = time_ns
    
    # Peptide N-terminus to receptor (closest CA) — use traj_ca coordinates
    if pep_n and receptor_ca:
        # Load pep_n and a sample of receptor CA into same traj
        dist_atoms2 = list(set(pep_n + receptor_ca))
        traj_dist2 = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=dist_atoms2, stride=stride)
        d2_map = {orig: i for i, orig in enumerate(dist_atoms2)}
        pn_idx = d2_map[pep_n[0]]
        rc_indices = [d2_map[i] for i in receptor_ca if i in d2_map]
        min_dists = []
        for frame in range(min(len(traj_dist2), len(time_ns))):
            pn_pos = traj_dist2.xyz[frame, pn_idx, :]
            ca_pos = traj_dist2.xyz[frame, rc_indices, :]
            dists = np.linalg.norm(ca_pos - pn_pos, axis=1)
            min_dists.append(np.min(dists))
        min_dists = np.array(min_dists)
    else:
        min_dists = None
    
    fig, axes = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
    if dist_bond is not None:
        axes[0].plot(time_ns_dist, dist_bond[:,0] * 10, color='tab:green', linewidth=0.8)
        axes[0].axhline(y=2.5, color='gray', linestyle='--', alpha=0.5, label='Bond limit (2.5 Å)')
        axes[0].set_ylabel('LYS-LNK Distance (Å)')
        axes[0].set_title(f'{name}: Key Structural Metrics')
        axes[0].legend()
    if min_dists is not None:
        time_ns_md = time_ns[:len(min_dists)]
        axes[1].plot(time_ns_md, min_dists * 10, color='tab:purple', linewidth=0.8)
        axes[1].set_ylabel('Peptide N-term to Receptor (Å)')
        axes[1].set_xlabel('Time (ns)')
    plt.tight_layout()
    plt.savefig(f"{OUT}/{name.lower().replace(' ', '_')}_distances.png")
    plt.close()
    
    if dist_bond is not None:
        np.savetxt(f"{OUT}/{name.lower().replace(' ', '_')}_bond_dist.dat", np.column_stack([time_ns_dist, dist_bond[:,0]*10]), header='Time(ns) Distance(A)')
    
    # === 4. Radius of gyration ===
    print("[4/8] Computing Rg...")
    rg_peptide = md.compute_rg(traj_ca.atom_slice(pc_idx))
    
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(time_ns, rg_peptide * 10, color='tab:red', linewidth=0.8)
    ax.axhline(y=np.mean(rg_peptide)*10, color='tab:red', linestyle='--', alpha=0.5)
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Peptide Rg (Å)')
    ax.set_title(f'{name}: Peptide Compactness')
    plt.tight_layout()
    plt.savefig(f"{OUT}/{name.lower().replace(' ', '_')}_rg.png")
    plt.close()
    
    # === 5. DSSP secondary structure ===
    print("[5/8] Computing DSSP...")
    try:
        traj_prot_all = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=receptor_ca + peptide_ca, stride=max(stride, total_frames//200))
        dssp = md.compute_dssp(traj_prot_all, simplified=True)
        
        # Count secondary structure types over time for peptide
        pep_dssp = dssp[:, len(receptor_ca):]
        ss_types = ['H', 'E', 'C']
        ss_counts = {ss: np.sum(pep_dssp == ss, axis=1) for ss in ss_types}
        
        fig, ax = plt.subplots(figsize=(8, 3))
        time_ss = np.arange(traj_prot_all.n_frames) * 10.0 * max(stride, total_frames//200) / 1000.0
        ax.stackplot(time_ss, ss_counts['H'], ss_counts['E'], ss_counts['C'],
                     labels=['Helix', 'Sheet', 'Coil'], colors=['#e74c3c', '#3498db', '#95a5a6'], alpha=0.8)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Residue Count')
        ax.set_title(f'{name}: Peptide Secondary Structure')
        ax.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig(f"{OUT}/{name.lower().replace(' ', '_')}_dssp.png")
        plt.close()
    except Exception as e:
        print(f"    DSSP skipped: {e}")
    
    # === 6. Membrane-specific analysis ===
    if is_membrane:
        print("[6/8] Membrane analysis...")
        
        # Load lipid P atoms
        p_atoms = [a.index for a in top.atoms if a.name == 'P' and a.residue.index >= 1140]
        if p_atoms:
            traj_p = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=p_atoms, stride=max(stride, total_frames//200))
            
            z0 = traj_p.xyz[0, :, 2]
            mid_z = np.median(z0)
            upper_idx = np.where(z0 > mid_z)[0]
            lower_idx = np.where(z0 < mid_z)[0]
            
            thickness = []
            for frame in range(traj_p.n_frames):
                z = traj_p.xyz[frame, :, 2]
                thickness.append(abs(np.mean(z[upper_idx]) - np.mean(z[lower_idx])))
            thickness = np.array(thickness)
            time_mem = np.arange(traj_p.n_frames) * 10.0 * max(stride, total_frames//200) / 1000.0
            
            fig, ax = plt.subplots(figsize=(8, 3))
            ax.plot(time_mem, thickness * 10, color='tab:green', linewidth=0.8)
            ax.axhline(y=np.mean(thickness)*10, color='tab:green', linestyle='--', alpha=0.5,
                       label=f'Mean: {np.mean(thickness)*10:.1f} Å')
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel('Bilayer Thickness (Å)')
            ax.set_title(f'{name}: Membrane Thickness (P-P)')
            ax.legend()
            plt.tight_layout()
            plt.savefig(f"{OUT}/{name.lower().replace(' ', '_')}_membrane_thickness.png")
            plt.close()
            
            np.savetxt(f"{OUT}/{name.lower().replace(' ', '_')}_thickness.dat", np.column_stack([time_mem, thickness*10]), header='Time(ns) Thickness(A)')
        
        # LNK tail z-position
        tail_atoms = [a.index for a in top.atoms if a.residue.index == 1139 and a.name in ['C55', 'C58', 'O56', 'O57', 'O59']]
        if tail_atoms:
            traj_tail = md.load_dcd(dcd_path, top=prmtop_path, atom_indices=tail_atoms, stride=max(stride, total_frames//200))
            tail_z = traj_tail.xyz[:, :, 2]
            time_tail = np.arange(traj_tail.n_frames) * 10.0 * max(stride, total_frames//200) / 1000.0
            
            fig, ax = plt.subplots(figsize=(8, 3))
            ax.plot(time_tail, tail_z * 10, color='tab:orange', linewidth=0.8, alpha=0.7)
            ax.axhline(y=0, color='gray', linestyle='--', alpha=0.3, label='Membrane center')
            ax.set_xlabel('Time (ns)')
            ax.set_ylabel('Tail z-position (Å)')
            ax.set_title(f'{name}: LNK Tail Position (diacid end)')
            ax.legend()
            plt.tight_layout()
            plt.savefig(f"{OUT}/{name.lower().replace(' ', '_')}_tail_z.png")
            plt.close()
    
    # === 7. Energy from log file ===
    print("[7/8] Plotting energies from log...")
    log_path = dcd_path.replace('.dcd', '.log').replace('membrane_rep1/', 'membrane_rep1/prod_v2.log')
    if 'membrane' in name.lower():
        log_path = f"{REPO}/exps/exp-F/md/membrane_rep1/prod_v2.log"
    else:
        log_path = f"{REPO}/exps/exp-F/md/rep1/prod.log"
    
    if os.path.exists(log_path):
        try:
            data = np.loadtxt(log_path, delimiter=',', skiprows=1)
            if data.ndim == 2 and data.shape[1] >= 9:
                step, t, pe, ke, te, temp, vol, dens, speed = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8]
                fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
                axes[0].plot(t/1000, temp, color='tab:red', linewidth=0.5)
                axes[0].set_ylabel('T (K)')
                axes[0].set_title(f'{name}: Thermodynamic Properties')
                axes[1].plot(t/1000, pe/1000, color='tab:blue', linewidth=0.5)
                axes[1].set_ylabel('PE (MJ/mol)')
                axes[2].plot(t/1000, dens, color='tab:green', linewidth=0.5)
                axes[2].set_ylabel('Density (g/mL)')
                axes[2].set_xlabel('Time (ns)')
                plt.tight_layout()
                plt.savefig(f"{OUT}/{name.lower().replace(' ', '_')}_energy.png")
                plt.close()
        except Exception as e:
            print(f"    Energy plot skipped: {e}")
    
    # === 8. Summary statistics ===
    print("[8/8] Writing summary...")
    with open(f"{OUT}/{name.lower().replace(' ', '_')}_summary.txt", 'w') as f:
        f.write(f"Analysis Summary: {name}\n")
        f.write(f"Trajectory length: {time_ns[-1]:.1f} ns ({traj_ca.n_frames} frames)\n")
        f.write(f"\nReceptor CA RMSD:\n")
        f.write(f"  Mean: {np.mean(rmsd_receptor)*10:.2f} Å\n")
        f.write(f"  Std:  {np.std(rmsd_receptor)*10:.2f} Å\n")
        f.write(f"  Max:  {np.max(rmsd_receptor)*10:.2f} Å\n")
        f.write(f"\nPeptide CA RMSD:\n")
        f.write(f"  Mean: {np.mean(rmsd_peptide)*10:.2f} Å\n")
        f.write(f"  Std:  {np.std(rmsd_peptide)*10:.2f} Å\n")
        f.write(f"  Max:  {np.max(rmsd_peptide)*10:.2f} Å\n")
        if dist_bond is not None:
            f.write(f"\nLYS-LNK Amide Bond:\n")
            f.write(f"  Mean: {np.mean(dist_bond)*10:.2f} Å\n")
            f.write(f"  Max:  {np.max(dist_bond)*10:.2f} Å\n")
            f.write(f"  Frames > 2.5 Å: {np.sum(dist_bond > 0.25)}/{len(dist_bond)} ({100*np.sum(dist_bond > 0.25)/len(dist_bond):.1f}%)\n")
        if min_dists is not None:
            f.write(f"\nPeptide N-term to Receptor:\n")
            f.write(f"  Mean: {np.mean(min_dists)*10:.2f} Å\n")
            f.write(f"  Initial: {min_dists[0]*10:.2f} Å\n")
            f.write(f"  Final: {min_dists[-1]*10:.2f} Å\n")
        f.write(f"\nPeptide Radius of Gyration:\n")
        f.write(f"  Mean: {np.mean(rg_peptide)*10:.2f} Å\n")
        if is_membrane and p_atoms:
            f.write(f"\nMembrane Thickness:\n")
            f.write(f"  Mean: {np.mean(thickness)*10:.1f} Å\n")
            f.write(f"  Std:  {np.std(thickness)*10:.1f} Å\n")
    
    elapsed = time.time() - t0
    print(f"Analysis complete in {elapsed/60:.1f} min")
    print(f"Output: {OUT}/{name.lower().replace(' ', '_')}*.png/dat")

# Run analyses
analyze("Solvent System", f"{REPO}/exps/exp-F/tleap/system.prmtop",
        f"{REPO}/exps/exp-F/md/rep1/prod.dcd", is_membrane=False)

analyze("Membrane System", f"{REPO}/exps/exp-F/membrane_build/system_final.prmtop",
        f"{REPO}/exps/exp-F/md/membrane_rep1/prod_v2.dcd", is_membrane=True)

print("\n" + "="*60)
print("ALL ANALYSES COMPLETE")
print(f"Output directory: {OUT}")
print("="*60)
