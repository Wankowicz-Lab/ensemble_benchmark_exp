# Bioemu Inference Runner
# This script will run Bioemu Inference

# Arguments: {pdb_id} {structure_fasta | structure} {num_samples}
# Outputs a merged/multi-conformer PDB file to provided output path


import sys
import os
import pandas as pd
import shutil

from bioemu.sample import main as sample
from scripts.helpers import dcd_to_pdb
import mdtraj as md


# merges the xtc samples and topology into a single PDB file with multiple conformations
def makeABigPDBFile(pdbId, tmp_path, output_path):
    xtc_file = os.path.join(tmp_path, 'samples.xtc')
    topology_file = os.path.join(tmp_path, 'topology.pdb')
    
    if not os.path.exists(xtc_file) or not os.path.exists(topology_file):
        print(f"[bioemu: inference.py] Files not found for {pdbId} in {tmp_path}. Skipping.")
        return None

    dcd_to_pdb.dcd_to_pdb(xtc_file, topology_file, output_path)
    print(f"[bioemu: inference.py] Saved combined PDB for {pdbId} at {output_path}")
    return output_path


def runBioemu(pdb_id, output_path, structure_fasta, num_samples):
    print(f"[bioemu: inference.py] Processing {pdb_id}...")
    
    temp_output_dir = os.path.join('./temp', pdb_id.lower())
    os.makedirs(temp_output_dir, exist_ok=True)
    
    sample(sequence=structure_fasta, 
           num_samples=num_samples, 
           output_dir=temp_output_dir)

    makeABigPDBFile(pdb_id, temp_output_dir, output_path)
    shutil.copyfile(os.path.join(temp_output_dir, 'samples.xtc'), './PDBs/' + pdb_id.lower() + '/bioemu_bin/bioemu.xtc')
    shutil.rmtree(temp_output_dir, ignore_errors=True)

if __name__ == "__main__":
    print("[bioemu: inference.py] Starting inference with bioemu")
    
    if len(sys.argv) < 5:
        print("[bioemu: inference.py] Usage: python inference.py <pdb_id> <output_path> <structure_fasta> <num_samples>")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    output_path = sys.argv[2]
    structure_fasta = sys.argv[3]
    num_samples = int(sys.argv[4])
    
    runBioemu(pdb_id, output_path, structure_fasta, num_samples)
    print(f"[bioemu: inference.py] Inference completed. Output saved to {output_path}")