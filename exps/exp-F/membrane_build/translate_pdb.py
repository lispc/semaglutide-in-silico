#!/usr/bin/env python3
import sys

def translate_pdb(input_path, output_path, dx=0, dy=0, dz=0):
    with open(input_path) as f:
        lines = f.readlines()
    
    with open(output_path, 'w') as f:
        for line in lines:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                x = float(line[30:38]) + dx
                y = float(line[38:46]) + dy
                z = float(line[46:54]) + dz
                f.write(f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}")
            else:
                f.write(line)

if __name__ == '__main__':
    translate_pdb(sys.argv[1], sys.argv[2], 0, 0, float(sys.argv[3]))
