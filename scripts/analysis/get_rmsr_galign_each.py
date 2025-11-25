# get_rmsr_galign_each_each.py <pdb_id>
# adds a rmsr_galign_each.csv to {PDB}/analysis as predictor|model|residue|residue_aa|RMSR

import numpy as np
import argparse
import sys
import os
import pandas as pd
from Bio.PDB import PDBParser
import pymol
from pymol import cmd

pdbparse = PDBParser(QUIET=True)

def align_with_pymol(reference_pdb, mobile_pdb, output_pdb):
    pymol.finish_launching(['pymol', '-qc'])
    
    
    import gc
    gc.collect()

    cmd.load(reference_pdb, "reference")
    cmd.load(mobile_pdb, "mobile")
    
    align_result = cmd.align("mobile", "reference")
    alignment_rmsd = align_result[0]
    
    cmd.save(output_pdb, "mobile", state=0)
    cmd.remove("all")
    cmd.delete("all")

    cmd.reinitialize()
    gc.collect()

    import time 
    time.sleep(1)
    
    
    print(f"[get_rmsr_galign_each.py] Aligned {mobile_pdb} to {reference_pdb} with RMSD {alignment_rmsd:.3f}")
    return alignment_rmsd

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


def get_rmsr_galign_each(pdb_id):
    PDB_FOLDER = f"./PDBs/{pdb_id.lower()}"
    deposited = f"{PDB_FOLDER}/{pdb_id.lower()}_final.pdb"

    structure = pdbparse.get_structure(pdb_id, deposited)
    multiconformer_centroids = {}
    for model in structure:
        if len(model) == 0:
            print(f"[get_rmsr_galign_each.py] Warning: No chains found in {deposited}, skipping...")
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

    print(f"[get_rmsr_galign_each.py] Collected {len(multiconformer_centroids)} multiconformer residue centroids from {deposited}.")

    predictors = ['bioemu', 'alphaflow', 'sam2', 'boltz2', 'openfold']
    RMSR_values = pd.DataFrame(columns=['predictor', 'model', 'residue', 'residue_aa', 'RMSR'])

    for i, predictor in enumerate(predictors):
        ensemble_path = f"{PDB_FOLDER}/{pdb_id.lower()}_{predictor}.pdb" # the ensemble path
        aligned_ensemble_path = f"{PDB_FOLDER}/{predictor}_bin/{pdb_id.lower()}_ensemble_pymol_aligned.pdb"
        
        print(f"[get_rmsr_galign_each.py] Processing {ensemble_path}...")

        if not os.path.exists(ensemble_path):
            print(f"[get_rmsr_galign_each.py] Warning: {ensemble_path} not found, skipping...")
            continue
            
        try:
            align_rmsd = align_with_pymol(deposited, ensemble_path, aligned_ensemble_path)
            print(f"[get_rmsr_galign_each.py] Successfully aligned {predictor} ensemble with RMSD {align_rmsd:.3f}")
        except Exception as e:
            print(f"[get_rmsr_galign_each.py] Error aligning {predictor} ensemble: {str(e)}")
            print(f"[get_rmsr_galign_each.py] Falling back to original ensemble")
            aligned_ensemble_path = ensemble_path

        ensemble = pdbparse.get_structure(f"{pdb_id}_{predictor}", aligned_ensemble_path)
      
        first_model = next(ensemble.get_models())
        if len(first_model) == 0:
            print(f"[get_rmsr_galign_each.py] Warning: No chains found in {aligned_ensemble_path}, skipping...")
            continue
            
        first_chain = next(first_model.get_chains())
        
        residue_names = {}
        res_ind = 0
        for residue in first_chain:
            residue_names[res_ind] = residue.resname
            res_ind += 1
            
        shift_by = 0
        max_shift_to_try = 5  # Try up to 5 positions shift
        best_shift = 0
        best_match_count = 0
        
        for test_shift in range(max_shift_to_try):
            match_count = 0
            for res_idx in range(min(10, len(residue_names))):  # Test with first 10 residues or fewer
                if res_idx >= len(residue_names):
                    break
                test_id = (res_idx + test_shift, residue_names[res_idx])
                if test_id in multiconformer_centroids:
                    match_count += 1
            
            if match_count > best_match_count:
                best_match_count = match_count
                best_shift = test_shift
                
        shift_by = best_shift
        print(f"[get_rmsr_galign_each.py] Using initial shift of {shift_by} for {predictor} ensemble")
        
        for model_idx, model in enumerate(ensemble):
            chain = next(model.get_chains())
            res_ind = 0
            current_shift = shift_by  # Each model starts with the base shift
            consecutive_mismatches = 0
            
            while res_ind < len(residue_names):
                residue = list(chain)[res_ind] if res_ind < len(list(chain)) else None
                if residue is None:
                    break
                    
                centroid = get_residue_centroid(residue)
                if centroid is None:
                    res_ind += 1
                    continue
                    
                residue_name = residue_names[res_ind]
                residue_id = (res_ind + current_shift, residue_name)
                
                # If this residue doesn't match, try additional shifts to find a match
                if residue_id not in multiconformer_centroids:
                    # Try looking ahead to see if there's a deletion
                    found_match = False
                    for test_shift in range(1, 4):  # Try up to 3 additional shifts
                        test_id = (res_ind + current_shift + test_shift, residue_name)
                        if test_id in multiconformer_centroids:
                            # Found a match with additional shift - adjust current_shift
                            print(f"[get_rmsr_galign_each.py] Detected possible missing residue at position {res_ind}, adjusting shift from {current_shift} to {current_shift + test_shift}")
                            current_shift += test_shift
                            found_match = True
                            consecutive_mismatches = 0
                            break
                            
                    if not found_match:
                        consecutive_mismatches += 1
                        if consecutive_mismatches > 5:
                            print(f"[get_rmsr_galign_each.py] Too many consecutive mismatches, skipping remaining residues in model {model_idx+1}")
                            break
                        res_ind += 1
                        continue
                else:
                    consecutive_mismatches = 0
                
                residue_id = (res_ind + current_shift, residue_name)
                if residue_id in multiconformer_centroids:
                    multiconformer_centroid = multiconformer_centroids[residue_id]
                    
                    diff = centroid - multiconformer_centroid
                    dist = np.sqrt(np.sum(diff**2))
                    
                    RMSR_values.loc[len(RMSR_values)] = [
                        predictor, 
                        model_idx + 1,  # Model number (1-indexed)
                        res_ind + current_shift, 
                        residue_name, 
                        dist
                    ]
                
                res_ind += 1

    output_file = f"{PDB_FOLDER}/analysis/rmsr_galign_each.csv"
    RMSR_values.to_csv(output_file, index=False)
    print(f"[get_rmsr_galign_each.py] Per-model RMSR values saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RMSR values for residues in predicted ensembles")
    parser.add_argument("pdb_id", help="PDB ID to process")

    args = parser.parse_args()
    pdb_id = args.pdb_id
    
    if not pdb_id.isalnum() or len(pdb_id) != 4:
        print("Error: PDB ID is wrong >> " + pdb_id)
        sys.exit(1)

    get_rmsr_galign_each(pdb_id)