#!/usr/bin/env python3
"""Full semaglutide MD: HSA + GLP-1 peptide + linker-FA with NZ-C8 restraint."""
import sys, os, argparse, time, numpy as np
import openmm as mm, openmm.app as app, openmm.unit as unit
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_C = f"{REPO}/exps/exp-C"

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--replica", type=int, default=1)
    p.add_argument("--gpu", type=str, default="0")
    p.add_argument("--nsteps", type=int, default=50_000_000)
    return p.parse_args()

def run(replica=1, gpu="0", nsteps=50_000_000):
    md_dir = f"{EXP_C}/md/full_sema/rep{replica}"
    os.makedirs(md_dir, exist_ok=True)

    prmtop = f"{EXP_C}/tleap/full_sema_hsa.prmtop"
    inpcrd = f"{EXP_C}/tleap/full_sema_hsa.inpcrd"
    system_xml = f"{EXP_C}/tleap/full_sema_hsa_system.xml"

    print(f"Loading {prmtop}")
    amber = pmd.load_file(prmtop, inpcrd)
    # Use get_coordinates() for plain numpy array (avoids Quantity issues)
    raw_positions = amber.get_coordinates()[0]  # shape (N, 3), units: angstroms

    # Find Lys26 NZ and LFA C8 in combined topology
    # HSA residues come first, then peptide, then LFA
    # Find the Lys in the peptide portion (residue name LYS after HSA residues)
    nz_idx = None; c8_idx = None
    for atom in amber.atoms:
        if atom.residue.name == 'LFA' and atom.name == 'C8':
            c8_idx = atom.idx
    # Find Lys26: LYS in the peptide range (pep_res_start to lfa_res_idx-1)
    # Lys26 is residue 19 in the 31-residue peptide (0-based), so pep_res_start + 19
    lfa_res_idx = None
    for r in amber.residues:
        if r.name == 'LFA': lfa_res_idx = r.idx; break
    pep_res_start = lfa_res_idx - 31
    lys26_res_idx = pep_res_start + 19  # Lys26 is peptide residue 19 (0-based)
    for atom in amber.atoms:
        if atom.residue.idx == lys26_res_idx and atom.name == 'NZ':
            nz_idx = atom.idx
            print(f"Found Lys26 NZ: res idx={atom.residue.idx}, atom idx={atom.idx}")
            break

    if nz_idx is None or c8_idx is None:
        print(f"ERROR: NZ={nz_idx}, C8={c8_idx}")
        return

    nz_pos = raw_positions[nz_idx]
    c8_pos = raw_positions[c8_idx]
    dist = np.linalg.norm(nz_pos - c8_pos) * 10
    print(f"NZ-C8 distance: {dist:.1f} A")

    # Position peptide: translate so Lys26 NZ is near LFA C8
    # Find peptide: residues immediately before LFA in the combined system
    pep_res_start = lfa_res_idx - 31  # GLP-1 has 31 residues
    pep_atom_indices = [a.idx for a in amber.atoms
                        if a.residue.idx >= pep_res_start and a.residue.idx < lfa_res_idx]
    print(f"Peptide: {len(pep_atom_indices)} atoms (residues {pep_res_start}-{lfa_res_idx-1})")

    # HSA atoms are everything before the peptide (non-water)
    hsa_last_idx = pep_res_start - 1
    hsa_atoms = [a.idx for a in amber.atoms
                 if a.residue.idx <= hsa_last_idx and a.residue.name not in ('WAT','HOH','SOL','Na+','Cl-')]

    # Translation: move NZ to within 1.5 A of C8
    translation = c8_pos - nz_pos + np.array([0.15, 0.0, 0.0])
    for idx in pep_atom_indices:
        raw_positions[idx] = raw_positions[idx] + translation

    # Verify
    new_nz = raw_positions[nz_idx]
    new_dist = np.linalg.norm(new_nz - c8_pos) * 10
    print(f"After translation: NZ-C8 distance = {new_dist:.1f} A")

    # Check min peptide-HSA distance
    min_d = float('inf')
    for pi in pep_atom_indices[::5]:
        p = raw_positions[pi]
        for hi in hsa_atoms[::20]:
            d = np.linalg.norm(p - raw_positions[hi])
            if d < min_d: min_d = d
    print(f"Min peptide-HSA distance: {min_d*10:.1f} A")

    # Build system fresh (caching with bond force causes ownership issues)
    print("Building system...")
    system = amber.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometers,
                                 constraints=app.HBonds, rigidWater=True)

    # Add CA restraint on HSA
    ca_force = mm.CustomExternalForce("5.0 * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    ca_force.addPerParticleParameter("x0"); ca_force.addPerParticleParameter("y0"); ca_force.addPerParticleParameter("z0")
    ca_count = 0
    for atom in amber.atoms:
        if atom.name == 'CA' and atom.residue.idx < pep_res_start:
            xyz = raw_positions[atom.idx].tolist()
            ca_force.addParticle(atom.idx, xyz)
            ca_count += 1
    system.addForce(ca_force)
    print(f"HSA CA restraint: {ca_count} atoms")

    # Add NZ-C8 soft harmonic restraint (k=500, gentle anchor)
    bond_force = mm.CustomBondForce("500.0 * (r - 0.145)^2")
    bond_force.addBond(nz_idx, c8_idx, [])
    system.addForce(bond_force)
    print(f"NZ-C8 restraint: {nz_idx}-{c8_idx} at 1.45 A, k=500 (soft anchor)")

    # Convert positions to nanometers for OpenMM
    positions_nm = raw_positions / 10.0  # A -> nm

    # Setup simulation
    integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
    integrator.setRandomNumberSeed(replica * 42)
    platform = mm.Platform.getPlatformByName('CUDA')
    simulation = app.Simulation(amber.topology, system, integrator, platform,
                                {'CudaDeviceIndex': gpu, 'CudaPrecision': 'mixed'})

    simulation.reporters.append(app.DCDReporter(f"{md_dir}/full_sema_traj.dcd", 50000))
    simulation.reporters.append(app.StateDataReporter(f"{md_dir}/full_sema_log.txt", 10000,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True, temperature=True, speed=True))
    simulation.reporters.append(app.CheckpointReporter(f"{md_dir}/full_sema_checkpoint.chk", 500000))

    simulation.context.setPositions(positions_nm)

    # Gentle NVT-only heating with position restraints, then add barostat
    heat_rf = mm.CustomExternalForce("200.0 * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    heat_rf.addPerParticleParameter("x0"); heat_rf.addPerParticleParameter("y0"); heat_rf.addPerParticleParameter("z0")
    for atom in amber.atoms:
        if atom.residue.idx <= hsa_last_idx or (atom.residue.idx >= pep_res_start and atom.residue.idx < lfa_res_idx):
            heat_rf.addParticle(atom.idx, raw_positions[atom.idx].tolist())
    hr_idx = system.addForce(heat_rf)
    print("Gentle heating NVT: 0→100K (50ps)")
    integrator.setTemperature(100 * unit.kelvin); simulation.step(25000)
    print("Gentle heating NVT: 100→200K (50ps)")
    integrator.setTemperature(200 * unit.kelvin); simulation.step(25000)
    print("Gentle heating NVT: 200→310K (50ps)")
    integrator.setTemperature(310 * unit.kelvin); simulation.step(25000)
    system.removeForce(hr_idx)
    print("NPT eq with barostat (100 ps, 310 K)")
    system.addForce(mm.MonteCarloBarostat(1*unit.bar, 310*unit.kelvin))
    simulation.step(50000)
    print(f"Production: {nsteps} steps ({nsteps*2e-6:.0f} ns)")
    sim_start = time.time(); steps_done = 0
    while steps_done < nsteps:
        chunk = min(500000, nsteps - steps_done)
        simulation.step(chunk); steps_done += chunk
        elapsed = time.time() - sim_start
        ns_done = steps_done * 2e-6; ns_day = ns_done / (elapsed / 86400)
        remaining = (nsteps - steps_done) * 2e-6 / ns_day * 24 if ns_day > 0 else 0
        print(f"[{time.strftime('%H:%M:%S')}] {ns_done:.0f}/{nsteps*2e-6:.0f} ns ({ns_day:.0f} ns/d, ~{remaining:.0f}h)")

    simulation.saveCheckpoint(f"{md_dir}/full_sema_final.chk")
    print(f"Done! {((time.time()-sim_start)/3600):.1f}h")

if __name__ == "__main__":
    args = parse_args()
    run(args.replica, args.gpu, args.nsteps)
