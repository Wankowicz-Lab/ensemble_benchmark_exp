# make_static_pdb.py <PDB_ID>
# For use with ./scripts/models/complete_alignment.sh

# MDTraj was failing to align because of atoms with non 1 occupancy.
# It was interpreting these as two different atoms.

# This script will remove atom lines with occupancy less than 0.5.

# Output will be from {pdb}_nowat.pdb to {pdb}_nowat_static.pdb 

import os
import sys

pdb_id = sys.argv[1].lower()
pdb_dir = "./PDBs/" + pdb_id;
nowat_pdb = pdb_dir + "/" + pdb_id + "_nowat.pdb"
static_pdb = pdb_dir + "/" + pdb_id + "_nowat_static.pdb"

if not os.path.exists(nowat_pdb):
    print(f"Error: Input file {nowat_pdb} not found.")
    sys.exit(1)

with open(nowat_pdb, 'r') as infile, open(static_pdb, 'w') as outfile:
    for line in infile:
        if line.startswith(('ATOM', 'HETATM')):
            try:
                occupancy = float(line[54:60].strip())
                if occupancy >= 0.5:
                    outfile.write(line)
            except (ValueError, IndexError):
                outfile.write(line)
        else:
            outfile.write(line)

print(f"[make_static_pdb.py] Processed {nowat_pdb} -> {static_pdb} by removing low occupancy atoms.")