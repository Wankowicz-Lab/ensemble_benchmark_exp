# get_rmsf.py <pdb_id>
# adds a rmsf.csv to {PDB}/analysis as predictor|residue|residue_aa|rmsf
import numpy as np
import argparse
import sys
import os
import mdtraj as md
import pandas as pd


# Stephanie Wankowicz: ensemble_bioinformatic_toolkit
def calculate_rmsf_from_residue_coords(residue_coords):
    mean_coords = np.mean(residue_coords, axis=0)
    rmsf = np.sqrt(np.mean(np.sum((residue_coords - mean_coords) ** 2, axis=1)))
    return rmsf

def get_rmsf_pdb(pdb_id):
    PDB_FOLDER = f"./PDBs/{pdb_id}"
    ANALYSIS_PATH = f"{PDB_FOLDER}/analysis"
    predictors = ['bioemu', 'alphaflow', 'sam2', 'boltz2', 'openfold']
    rmsf_values = []

    for i, predictor in enumerate(predictors):
        ensemble_path = f"{PDB_FOLDER}/{pdb_id.lower()}_{predictor}.pdb"
        print(f"[get_rmsf.py] Processing {ensemble_path}...")

        if not os.path.exists(ensemble_path):
            print(f"[get_rmsf.py] Warning: {ensemble_path} not found, skipping...")
            continue

        frames = md.load(ensemble_path)
        residues = frames.topology.residues
        
        for residue in residues:
            residue_indices = [atom.index for atom in residue.atoms]
            residue_coords = frames.xyz[:, residue_indices, :]
            rmsf = calculate_rmsf_from_residue_coords(residue_coords)
            rmsf_values.append((predictor, residue.index, residue.name, rmsf))
        print(f"[get_rmsf.py] Successfully calculated RMSF for residues in {predictor} ensemble.")

    new_df = pd.DataFrame(rmsf_values, columns=['predictor', 'residue', 'residue_aa', 'rmsf'])
    if not os.path.exists(ANALYSIS_PATH):
        os.makedirs(ANALYSIS_PATH, exist_ok=True)

    csv_path = f"{ANALYSIS_PATH}/rmsf.csv"
    new_df.to_csv(csv_path, index=False)
    print(f"[get_rmsf.py] RMSF values saved to {csv_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RMSF values for residues in ensembles made by predictors")
    parser.add_argument("pdb_id", help="PDB ID to process")

    args = parser.parse_args()
    pdb_id = args.pdb_id
    
    if not pdb_id.isalnum() or len(pdb_id) != 4:
        print("Error: PDB ID is wrong >> " + pdb_id)
        sys.exit(1)

    get_rmsf_pdb(pdb_id)

    