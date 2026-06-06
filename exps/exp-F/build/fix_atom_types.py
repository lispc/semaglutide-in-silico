#!/usr/bin/env python3
"""Fix GAFF2 atom types in LNK mol2."""
import sys

def fix_types(input_path, output_path):
    with open(input_path) as f:
        lines = f.readlines()
    
    atoms = {}
    bonds = []
    in_atom = False
    in_bond = False
    
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atom = True
            continue
        elif line.startswith("@<TRIPOS>BOND"):
            in_atom = False
            in_bond = True
            continue
        elif line.startswith("@<TRIPOS>") and in_bond:
            in_bond = False
        
        if in_atom:
            parts = line.split()
            if len(parts) >= 9:
                atoms[int(parts[0])] = {
                    'name': parts[1], 'type': parts[5],
                    'x': float(parts[2]), 'y': float(parts[3]), 'z': float(parts[4]),
                    'charge': float(parts[8])
                }
        elif in_bond:
            parts = line.split()
            if len(parts) >= 4:
                bonds.append((int(parts[1]), int(parts[2]), parts[3]))
    
    # Build neighbor map
    neighbors = {i: [] for i in atoms}
    for a1, a2, _ in bonds:
        neighbors[a1].append(a2)
        neighbors[a2].append(a1)
    
    # Determine fixes based on ORIGINAL types
    type_fixes = {}
    for idx, a in atoms.items():
        nbs = neighbors[idx]
        nb_types = [atoms[j]['type'] for j in nbs]
        
        # Carbonyl carbon: connected to O (carbonyl) and (N or C)
        if a['type'] == 'c3' and 'o' in nb_types:
            for j in nbs:
                if atoms[j]['type'] == 'o' and len(neighbors[j]) == 1:
                    type_fixes[idx] = 'c'
                    break
            continue
        
        # Ether oxygen: connected to two sp3 carbons
        if a['type'] == 'o' and len(nbs) == 2 and all(atoms[j]['type'] == 'c3' for j in nbs):
            type_fixes[idx] = 'os'
            continue
    
    # Apply fixes
    fixes = []
    for idx, new_type in type_fixes.items():
        old = atoms[idx]['type']
        atoms[idx]['type'] = new_type
        fixes.append(f"{idx} {atoms[idx]['name']}: {old} -> {new_type}")
    
    print(f"Fixed {len(fixes)} atom types:")
    for f in fixes:
        print(f"  {f}")
    
    # Write output
    with open(output_path, 'w') as f:
        f.write(f"@<TRIPOS>MOLECULE\nLNK\n  {len(atoms)}   {len(bonds)}     0     0     0\n")
        f.write("SMALL\nGAFF2\nlinker-C18-noN\n\n@<TRIPOS>ATOM\n")
        for idx in sorted(atoms):
            a = atoms[idx]
            f.write(f"{idx:>4} {a['name']:<6} {a['x']:>10.4f} {a['y']:>10.4f} {a['z']:>10.4f} {a['type']:>4}   1 LNK {a['charge']:>10.6f}\n")
        f.write("\n@<TRIPOS>BOND\n")
        for i, (a1, a2, btype) in enumerate(bonds, 1):
            f.write(f"{i:>4} {a1:>4} {a2:>4} {btype}\n")
        f.write("\n@<TRIPOS>SUBSTRUCTURE\n     1 LNK         1 TEMP              0 ****  ****    0 ROOT\n")
    
    print(f"Wrote {output_path}")

if __name__ == '__main__':
    fix_types(sys.argv[1], sys.argv[2])
