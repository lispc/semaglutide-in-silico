import math
import csv
import os

try:
    from scipy.spatial import cKDTree
    HAS_KDTREE = True
except:
    HAS_KDTREE = False

os.chdir('/home/scroll/personal/semaglutide-in-silico/exps/exp-F/membrane_build')

def read_mapping(csv_file):
    mapping = {}
    with open(csv_file) as f:
        f.readline()
        reader = csv.DictReader(f)
        for row in reader:
            res = row['residue'].strip()
            search = row['search'].strip()
            replace = row['replace'].strip()
            order = int(row['order'])
            ter = row['TER'].strip() == 'True'
            num_atom = int(row['num_atom'])
            if res not in mapping:
                mapping[res] = {'atoms': {}, 'ter': ter, 'num_atom': num_atom}
            mapping[res]['atoms'][search] = {'replace': replace, 'order': order}
    return mapping

mapping = read_mapping('/home/scroll/miniforge3/envs/cgas-md/lib/python3.11/site-packages/packmol_memgen/lib/charmmlipid2amber/charmmlipid2amber.csv')
print(f"Loaded mappings for: {list(mapping.keys())[:5]}...")

# Read protein
protein_atoms = []
with open('protein.pdb') as f:
    for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            protein_atoms.append({
                'record': 'ATOM',
                'atomname': line[12:16].strip(),
                'resname': line[17:20].strip(),
                'chain': line[21],
                'resnum': int(line[22:26].strip()),
                'x': float(line[30:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54]),
                'element': line[77:78].strip() if len(line) > 77 else ''
            })

protein_chains = {}
for atom in protein_atoms:
    c = atom['chain']
    if c not in protein_chains:
        protein_chains[c] = []
    protein_chains[c].append(atom)

print(f"Protein: {len(protein_atoms)} atoms, chains: {list(protein_chains.keys())}")

# Read membrane_only.pdb and group by TER
mem_molecules = []
current_mol = []
with open('membrane_only.pdb') as f:
    for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            current_mol.append({
                'record': 'ATOM',
                'atomname': line[12:16].strip(),
                'resname': line[17:20].strip(),
                'chain': line[21],
                'resnum': int(line[22:26].strip()),
                'x': float(line[30:38]),
                'y': float(line[38:46]),
                'z': float(line[46:54]),
                'element': line[77:78].strip() if len(line) > 77 else ''
            })
        elif line.startswith('TER') and current_mol:
            mem_molecules.append(current_mol)
            current_mol = []
    if current_mol:
        mem_molecules.append(current_mol)

print(f"Membrane molecules: {len(mem_molecules)}")

lipid_mols = [m for m in mem_molecules if m[0]['resname'] in ['POP','CHL']]
water_mols = [m for m in mem_molecules if m[0]['resname'] == 'TIP']
ion_mols = [m for m in mem_molecules if m[0]['resname'] in ['K+','Cl-']]
print(f"Lipids: {len(lipid_mols)}, Waters: {len(water_mols)}, Ions: {len(ion_mols)}")

# Distance check
prot_coords = [(a['x'], a['y'], a['z']) for a in protein_atoms]
if HAS_KDTREE:
    kdtree = cKDTree(prot_coords)
    print("Using KD-tree")
else:
    print("Warning: using brute force")

def min_dist_to_protein(mol):
    if HAS_KDTREE:
        mol_coords = [(a['x'], a['y'], a['z']) for a in mol]
        mind, _ = kdtree.query(mol_coords, k=1)
        return mind.min()
    else:
        mind = float('inf')
        for atom in mol:
            ax, ay, az = atom['x'], atom['y'], atom['z']
            for px, py, pz in prot_coords:
                d = (ax-px)**2 + (ay-py)**2 + (az-pz)**2
                if d < mind:
                    mind = d
                    if mind < 6.25:
                        return math.sqrt(mind)
        return math.sqrt(mind)

print("Checking lipid overlap...")
keep_lipids = []
remove_count = 0
for i, mol in enumerate(lipid_mols):
    d = min_dist_to_protein(mol)
    if d < 2.5:
        remove_count += 1
    else:
        keep_lipids.append(mol)
    if i % 50 == 0:
        print(f"  {i}/{len(lipid_mols)} checked, removed {remove_count}")

print(f"Kept {len(keep_lipids)} lipids, removed {remove_count}")

def format_atomnum(n):
    if n > 99999:
        return f"{n:>5X}"
    return f"{n:>5d}"

def format_resnum(n):
    if n > 9999:
        return f"{n:>4X}"
    return f"{n:>4d}"

# Write merged PDB
chain_list = list('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')
atom_counter = 1

with open('merged_system.pdb', 'w') as out:
    # Protein
    for chain_name in sorted(protein_chains.keys()):
        atoms = protein_chains[chain_name]
        for atom in atoms:
            out.write(f"ATOM  {format_atomnum(atom_counter)} {atom['atomname']:>4s} {atom['resname']:>3s}{chain_name:>1s}{format_resnum(atom['resnum'])}    {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}  1.00  0.00           {atom['element']:>1s}\n")
            atom_counter += 1
        out.write("TER\n")
    
    # Lipids
    chain_idx = 0
    for mol in keep_lipids:
        chain = chain_list[chain_idx % len(chain_list)]
        resname0 = mol[0]['resname']
        map_key = {'POP':'POPC','CHL':'CHL1'}.get(resname0, resname0)
        
        if map_key in mapping:
            molmap = mapping[map_key]
            converted = []
            for atom in mol:
                search = f"{atom['atomname']:>4s} {map_key:>3s}"
                if search in molmap['atoms']:
                    rep = molmap['atoms'][search]['replace']
                else:
                    search2 = f"{atom['atomname']} {map_key}"
                    if search2 in molmap['atoms']:
                        rep = molmap['atoms'][search2]['replace']
                    else:
                        continue
                new_atomname = rep[:4].strip()
                new_resname = rep[4:].strip()
                order = molmap['atoms'][search]['order'] if search in molmap['atoms'] else molmap['atoms'][search2]['order']
                converted.append({
                    'atomname': new_atomname,
                    'resname': new_resname,
                    'order': order,
                    'x': atom['x'], 'y': atom['y'], 'z': atom['z'],
                    'element': atom['element']
                })
            
            converted.sort(key=lambda a: a['order'])
            
            last_resname = None
            resnum = 0
            for atom in converted:
                if atom['resname'] != last_resname:
                    resnum += 1
                    last_resname = atom['resname']
                out.write(f"ATOM  {format_atomnum(atom_counter)} {atom['atomname']:>4s} {atom['resname']:>3s}{chain:>1s}{format_resnum(resnum)}    {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}  1.00  0.00           {atom['element']:>1s}\n")
                atom_counter += 1
            out.write("TER\n")
        else:
            resnum = 1
            for atom in mol:
                out.write(f"ATOM  {format_atomnum(atom_counter)} {atom['atomname']:>4s} {atom['resname']:>3s}{chain:>1s}{format_resnum(resnum)}    {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}  1.00  0.00           {atom['element']:>1s}\n")
                atom_counter += 1
            out.write("TER\n")
        
        chain_idx += 1
    
    # Waters
    wchain_idx = 0
    wresnum = 1
    for i, mol in enumerate(water_mols):
        if wresnum > 9999:
            wchain_idx += 1
            wresnum = 1
        chain = chain_list[wchain_idx % len(chain_list)]
        for atom in mol:
            out.write(f"ATOM  {format_atomnum(atom_counter)} {atom['atomname']:>4s} {'WAT':>3s}{chain:>1s}{format_resnum(wresnum)}    {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}  1.00  0.00           {atom['element']:>1s}\n")
            atom_counter += 1
        wresnum += 1
    
    # Ions
    ichain_idx = 0
    iresnum = 1
    for i, mol in enumerate(ion_mols):
        if iresnum > 9999:
            ichain_idx += 1
            iresnum = 1
        chain = chain_list[ichain_idx % len(chain_list)]
        for atom in mol:
            out.write(f"ATOM  {format_atomnum(atom_counter)} {atom['atomname']:>4s} {atom['resname']:>3s}{chain:>1s}{format_resnum(iresnum)}    {atom['x']:>8.3f}{atom['y']:>8.3f}{atom['z']:>8.3f}  1.00  0.00           {atom['element']:>1s}\n")
            atom_counter += 1
        iresnum += 1
    
    out.write("END\n")

print(f"Wrote merged_system.pdb with {atom_counter-1} atoms")
