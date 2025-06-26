# get_secondary_structure.py <pdb_id>
# Uses DSSP to get secondary structure information from a PDB file.

# https://github.com/PDB-REDO/dssp

# INSTALL WITH:
# cd /home/aslamaj/workdir/ensemble_benchmark_2/bin/misc/dssp
# rm -rf build/
# mkdir build
# cd build

# # Configure with all paths pointing to local directory
# cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local \
#       -DCMAKE_INSTALL_LIBDIR=$HOME/.local/lib \
#       -DCMAKE_INSTALL_INCLUDEDIR=$HOME/.local/include \
#       -DCMAKE_INSTALL_BINDIR=$HOME/.local/bin \
#       -DCMAKE_INSTALL_DATADIR=$HOME/.local/share \
#       ..

# make
# cmake --install .



import argparse
import pandas as pd
import os
import subprocess
import tempfile
from Bio import PDB
from Bio.PDB.DSSP import DSSP

def get_secondary_structure_from_dssp(pdb_id, pdb_file_path):
    try:
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure(pdb_id, pdb_file_path)
        
        model = structure[0]
        
        # .local/bin/mkdssp
        dssp_exe = os.path.expanduser("~/.local/bin/mkdssp")

        os.environ['LIBCIFPP_DATA_DIR'] = os.path.realpath('./bin/libcifpp')

        dssp = DSSP(model, pdb_file_path, dssp=dssp_exe)
        
        ss_data = []
        for key in dssp.keys():
            chain_id = key[0]
            residue_id = key[1][1]
            residue_name = dssp[key][1]
            ss_structure = dssp[key][2]
            
            ss_mapping = {
                'H': 'H',  # Alpha helix
                'B': 'E',  # Beta bridge -> Beta sheet
                'E': 'E',  # Extended strand (beta sheet)
                'G': 'G',  # 3-10 helix
                'I': 'I',  # Pi helix
                'T': 'T',  # Turn
                'S': 'S',  # Bend
                ' ': 'C',  # Coil
                '-': 'C'   # Coil
            }
            
            ss_data.append({
                'residue': residue_id,
                'chain': chain_id,
                'amino_acid': residue_name,
                'secondary_structure': ss_mapping.get(ss_structure, 'C')
            })
        
        return pd.DataFrame(ss_data)
        
    except Exception as e:
        print(f"Error running DSSP: {e}")
        return None



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get secondary structure from PDB file.")
    parser.add_argument("pdb_id", type=str, help="PDB ID of the structure")
    args = parser.parse_args()

    pdb_id = args.pdb_id
    pdb_file_path = f"./PDBs/{pdb_id}/{pdb_id}_final.pdb"

    if not os.path.exists(pdb_file_path):
        print(f"PDB file {pdb_file_path} does not exist.")
        exit(1)

    ssdf = get_secondary_structure_from_dssp(pdb_id, pdb_file_path)

    # RESIDUE|SECONDARY_STRUCTURE
    newdf = ssdf[['residue', 'secondary_structure']].copy()
    newdf.to_csv(f"./PDBs/{pdb_id}/analysis/secondary_structure.csv", index=False)
    print(f"Secondary structure for {pdb_id} has been processed.")