#!/usr/bin/env python3
"""Load checkpoint and save coordinates without PBC wrapping."""
import openmm as mm
import openmm.app as app
import openmm.unit as unit
import parmed as pmd

REPO = "/home/scroll/personal/semaglutide-in-silico"
EXP_F = f"{REPO}/exps/exp-F"
PRMTOP = f"{EXP_F}/tleap/system.prmtop"
CHK = f"{EXP_F}/md/rep1/equil.chk"

def add_restraints(system, topology, positions):
    receptor_force = mm.CustomExternalForce("10.0 * periodicdistance(x, y, z, x0, y0, z0)^2")
    receptor_force.addPerParticleParameter("x0")
    receptor_force.addPerParticleParameter("y0")
    receptor_force.addPerParticleParameter("z0")

    pep_force = mm.CustomExternalForce("100.0 * periodicdistance(x, y, z, x0, y0, z0)^2")
    pep_force.addPerParticleParameter("x0")
    pep_force.addPerParticleParameter("y0")
    pep_force.addPerParticleParameter("z0")

    res_idx = 0
    for chain in topology.chains():
        for res in chain.residues():
            if res.name in ('WAT', 'Na+', 'Cl-', 'HOH', 'SOL', 'LNK'):
                res_idx += 1
                continue
            if 1113 <= res_idx <= 1118:
                for atom in res.atoms():
                    if atom.name in ('N', 'CA', 'C'):
                        xyz = positions[atom.index]
                        pep_force.addParticle(atom.index, [xyz.x, xyz.y, xyz.z])
                res_idx += 1
                continue
            if res_idx <= 1112:
                for atom in res.atoms():
                    if atom.name == 'CA':
                        xyz = positions[atom.index]
                        receptor_force.addParticle(atom.index, [xyz.x, xyz.y, xyz.z])
            res_idx += 1

    system.addForce(receptor_force)
    system.addForce(pep_force)
    return system

print("Loading system...")
amber = pmd.load_file(PRMTOP, f"{EXP_F}/tleap/system.inpcrd")
system = amber.createSystem(
    nonbondedMethod=app.PME,
    nonbondedCutoff=1.0*unit.nanometers,
    constraints=app.HBonds,
    rigidWater=True,
)

topology = app.AmberPrmtopFile(PRMTOP).topology
positions = amber.positions  # dummy, will be overwritten

system = add_restraints(system, topology, positions)

integrator = mm.LangevinIntegrator(310*unit.kelvin, 1.0/unit.picoseconds, 2.0*unit.femtoseconds)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}

system.addForce(mm.MonteCarloBarostat(1*unit.bar, 310*unit.kelvin))

simulation = app.Simulation(topology, system, integrator, platform, properties)

print(f"Loading checkpoint from {CHK}")
with open(CHK, 'rb') as f:
    simulation.context.loadCheckpoint(f.read())

print("Saving unwrapped coordinates...")
state = simulation.context.getState(getPositions=True, enforcePeriodicBox=False)
with open(f"{EXP_F}/md/rep1/equilibrated_unwrapped.pdb", "w") as f:
    app.PDBFile.writeFile(topology, state.getPositions(), f)

print("Saved equilibrated_unwrapped.pdb")
