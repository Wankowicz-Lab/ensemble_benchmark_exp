# get_bbx.py
# Get bounding box/sphere of deposited protein structure.
# Outputs X, Y, Z and Range of a sphere that fits around the structure.
# Does not need to be that exact, just enough that Phaser doesn't align to the wrong place.

import os
import sys
import numpy as np


def get_protein_bounding_box(pdb_path):
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")
    
    coords = []
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue
    
    if not coords:
        raise ValueError(f"No atom coordinates found in PDB file: {pdb_path}")
    
    coords = np.array(coords)
    
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    
    centroid = (min_coords + max_coords) / 2
    centroid_x, centroid_y, centroid_z = centroid
    
    ranges = max_coords - min_coords
    box_size = np.max(ranges) / 2
    
    return centroid_x, centroid_y, centroid_z, box_size


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python get_bbx.py <path_to_pdb_file>")
        sys.exit(1)
    
    pdb_file_path = sys.argv[1]
    
    try:
        x, y, z, box_size = get_protein_bounding_box(pdb_file_path)
        print(f"{x:.2f},{y:.2f},{z:.2f},{box_size:.2f}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)