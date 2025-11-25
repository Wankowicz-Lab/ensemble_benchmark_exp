# get_rmsr_mr.py <pdb_id>
# adds a rmsr_mr.csv to {PDB}/analysis as predictor|residue|residue_aa|RMSR

import numpy as np
import argparse
import sys
import os
import pandas as pd
from Bio.PDB import PDBParser
import pymol
from pymol import cmd

pdbparse = PDBParser(QUIET=True)

def get_multiconformer_residue_centroid(residue):
    altloc_groups = {}
    for atom in residue:
        altloc = atom.get_altloc()
        if altloc == ' ':
            altloc = 'A'
        
        occupancy = atom.get_occupancy()
        if occupancy is None or occupancy == 0.0:
            continue

        if atom.get_name() == 'H':
            continue

        if altloc not in altloc_groups:
            altloc_groups[altloc] = {'coords': [], 'occupancy': occupancy}
        
        altloc_groups[altloc]['coords'].append(atom.get_coord())

    if not altloc_groups:
        return None

    weighted_sum = np.zeros(3)
    total_weight = 0.0

    for altloc, data in altloc_groups.items():
        coords = np.array(data['coords'])
        occ = data['occupancy']
        centroid = coords.mean(axis=0)
        weighted_sum += occ * centroid
        total_weight += occ

    if total_weight == 0.0:
        return None

    return weighted_sum / total_weight


def get_residue_centroid(residue):
    coords = []
    for atom in residue:
        if atom.get_name() == 'H':
            continue
        coords.append(atom.get_coord())
    
    if not coords:
        return None
    
    return np.mean(coords, axis=0)


def get_rmsr_mr(pdb_id):
    PDB_FOLDER = f"./PDBs/{pdb_id.lower()}"
    deposited = f"{PDB_FOLDER}/{pdb_id.lower()}_final.pdb"

    structure = pdbparse.get_structure(pdb_id, deposited)
    multiconformer_centroids = {}
    for model in structure:
        if len(model) == 0:
            print(f"[get_rmsr_mr.py] Warning: No chains found in {deposited}, skipping...")
            continue
        chains = model.get_chains() 
        chain = next(chains) # for now, first chain only for consistency.
        res_ind = 0 # added bc residues are in order, but resnumbers may start at 0, 1, or 2.
        for residue in chain:
            centroid = get_multiconformer_residue_centroid(residue)
            if centroid is not None:
                residue_id = (res_ind, residue.resname)
                multiconformer_centroids[residue_id] = centroid
            res_ind += 1

    print(f"[get_rmsr_mr.py] Collected {len(multiconformer_centroids)} multiconformer residue centroids from {deposited}.")

    predictors = ['bioemu', 'alphaflow', 'sam2', 'boltz2', 'openfold']
    RMSR_values = pd.DataFrame(columns=['predictor', 'residue', 'residue_aa', 'RMSR'])

    for i, predictor in enumerate(predictors):
        ensemble_path = f"{PDB_FOLDER}/{pdb_id.lower()}_{predictor}.pdb" # the ensemble path
       
        print(f"[get_rmsr_mr.py] Processing {ensemble_path}...")

        if not os.path.exists(ensemble_path):
            print(f"[get_rmsr_mr.py] Warning: {ensemble_path} not found, skipping...")
            continue
            

        ensemble = pdbparse.get_structure(f"{pdb_id}_{predictor}", ensemble_path)
        
        first_model = next(ensemble.get_models())
        if len(first_model) == 0:
            print(f"[get_rmsr_mr.py] Warning: No chains found in {ensemble_path}, skipping...")
            continue
            
        first_chain = next(first_model.get_chains())
        
        residue_names = {}
        res_ind = 0
        for residue in first_chain:
            residue_names[res_ind] = residue.resname
            res_ind += 1
        
        residue_centroids = {}
        
        for model in ensemble:
            chain = next(model.get_chains())
            res_ind = 0
            for residue in chain:
                centroid = get_residue_centroid(residue)
                if centroid is not None:
                    if res_ind not in residue_centroids:
                        residue_centroids[res_ind] = []
                    residue_centroids[res_ind].append(centroid)
                res_ind += 1
        
        shift_by = 0
        for res_ind, centroids in residue_centroids.items():
            if not centroids:
                continue
                
            residue_name = residue_names[res_ind]
            residue_id = (res_ind + shift_by, residue_name)

            
            if residue_id not in multiconformer_centroids:
                print(f"[get_rmsr_mr.py] Warning: {residue_id} is not in the multiconformer")
                shift_by += 1
                continue
                
            multiconformer_centroid = multiconformer_centroids[residue_id]
            centroids_array = np.array(centroids)
            
            diffs = centroids_array - multiconformer_centroid
            RMSR = np.sqrt(np.mean(np.sum(diffs**2, axis=1)))
            
            RMSR_values.loc[len(RMSR_values)] = [
                predictor, 
                res_ind + shift_by, 
                residue_name, 
                RMSR
            ]

    RMSR_values.to_csv(f"{PDB_FOLDER}/analysis/rmsr_mr.csv", index=False)
    print(f"[get_rmsr_mr.py] RMSR values saved to {PDB_FOLDER}/analysis/rmsr_mr.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RMSR values for residues in predicted ensembles")
    parser.add_argument("pdb_id", help="PDB ID to process")

    args = parser.parse_args()
    pdb_id = args.pdb_id
    
    if not pdb_id.isalnum() or len(pdb_id) != 4:
        print("Error: PDB ID is wrong >> " + pdb_id)
        sys.exit(1)

    get_rmsr_mr(pdb_id)