# get_rmsf.py <pdb_id> --<mode, default=heavy | CA, backbone, heavy, all>
# adds a rmsf.csv to {PDB}/analysis as predictor|residue|residue_aa|rmsf
import numpy as np
import argparse
import sys
import os
import pandas as pd
from Bio.PDB import PDBParser, is_aa
from collections import defaultdict


def calculate_rmsf(residue_coords):
    mean_coords = np.mean(residue_coords, axis=0)
    rmsf = np.sqrt(np.mean(np.sum((residue_coords - mean_coords) ** 2, axis=1)))
    return rmsf

def atom_selector(atom, mode="heavy"):
    if mode == "CA":
        return atom.get_id() == "CA"
    if mode == "backbone":
        return atom.get_id() in {"N", "CA", "C", "O"}
    if mode == "heavy":
        return atom.element != "H"
    if mode == "all":
        return True
    return True


def get_rmsf_pdb(pdb_id, mode="heavy"):
    PDB_FOLDER = f"./PDBs/{pdb_id}"

    if (pdb_id.find('_')):
        pdb_id = pdb_id.split('_')[0]
    ANALYSIS_PATH = f"{PDB_FOLDER}/analysis"
    predictors = ['bioemu', 'alphaflow', 'sam2', 'boltz2', 'openfold']
    all_rmsf_dfs = [] 
    
    parser = PDBParser(QUIET=True)

    for predictor in predictors:
        ensemble_path = f"{PDB_FOLDER}/{pdb_id.lower()}_{predictor}.pdb"
        print(f"[get_rmsf.py] Processing {ensemble_path}...")

        if not os.path.exists(ensemble_path):
            print(f"[get_rmsf.py] Warning: {ensemble_path} not found, skipping...")
            continue

        structure = parser.get_structure(predictor, ensemble_path)
        
    

        # centroid calculation per residue:
        # residue_rmsf = {}
        # for model in structure:
        #     for chain in model:
        #         for residue in chain:
        #             if is_aa(residue):
        #                 res_id = (chain.id, residue.id[1])
        #                 if res_id not in residue_rmsf:
        #                     residue_rmsf[res_id] = []
        #                 atom_coords = [atom.coord for atom in residue if atom.element != 'H']
        #                 residue_rmsf[res_id].append(np.mean(atom_coords, axis=0))
        
        # rmsf_values = {res_id: calculate_rmsf(np.array(coords)) for res_id, coords in residue_rmsf.items()}
        


        # conventional calculation per atom:
        per_atom_traj = defaultdict(list)  # (chain_id, resseq, atom_name) -> list of 3D coords
        all_res_ids = set()
        n_frames = 0

        for model in structure:
            n_frames += 1
            for chain in model:
                for residue in chain:
                    if not is_aa(residue, standard=True):
                        continue
                    res_id = (chain.id, residue.id[1])
                    all_res_ids.add(res_id)
                    for atom in residue:
                        if atom_selector(atom, mode=mode):
                            key = (chain.id, residue.id[1], atom.get_id())
                            per_atom_traj[key].append(atom.coord)

        per_atom_rmsf = {}
        for key, coords in per_atom_traj.items():
            coords = np.asarray(coords)
            if coords.shape[0] == n_frames: # drop atoms that are not consistent / in every frame
                per_atom_rmsf[key] = calculate_rmsf(coords)
        
        residue_rmsf = defaultdict(list)
        for (chain_id, resseq, atom_name), rmsf in per_atom_rmsf.items():
            residue_rmsf[(chain_id, resseq)].append(rmsf)

        rmsf_values = {res_id: float(np.mean(vals)) for res_id, vals in residue_rmsf.items()}



        rmsf_df = pd.DataFrame(list(rmsf_values.items()), columns=['Residue', 'RMSF'])
        rmsf_df[['Chain', 'Resi']] = pd.DataFrame(rmsf_df['Residue'].tolist(), index=rmsf_df.index)
        rmsf_df = rmsf_df.drop(columns=['Residue'])
        rmsf_df['PDB'] = predictor
        
        residue_aa_dict = {}
        for res_id in rmsf_values.keys():
            chain_id, residue_idx = res_id
            for model in structure:
                for chain in model:
                    if chain.id == chain_id:
                        for residue in chain:
                            if is_aa(residue) and residue.id[1] == residue_idx:
                                residue_aa_dict[res_id] = residue.get_resname()
                                break
        
        rmsf_df['residue_aa'] = rmsf_df.apply(
            lambda row: residue_aa_dict.get((row['Chain'], row['Resi']), ''), axis=1
        )
        
        all_rmsf_dfs.append(rmsf_df)
        
        if not os.path.exists(ANALYSIS_PATH):
            os.makedirs(ANALYSIS_PATH, exist_ok=True)
        #rmsf_df.to_csv(f"{ANALYSIS_PATH}/rmsf_output_{predictor}.csv", index=False)
        
        print(f"[get_rmsf.py] Successfully calculated RMSF for residues in {predictor} ensemble.")
    
    if all_rmsf_dfs:
        concatenated_rmsf_df = pd.concat(all_rmsf_dfs, ignore_index=True)
        
        final_df = pd.DataFrame({
            'predictor': concatenated_rmsf_df['PDB'],
            'residue': concatenated_rmsf_df['Resi'],
            'residue_aa': concatenated_rmsf_df['residue_aa'],
            'rmsf': concatenated_rmsf_df['RMSF']
        })
        
        filename = "rmsf.csv"
        if (mode != "heavy"):
            filename = f"rmsf_{mode}.csv"
        csv_path = f"{ANALYSIS_PATH}/{filename}"
        final_df.to_csv(csv_path, index=False)
        print(f"[get_rmsf.py] RMSF values saved to {csv_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate RMSF values for residues in ensembles made by predictors")
    parser.add_argument("pdb_id", help="PDB ID to process")
    parser.add_argument("--mode", choices=["CA", "backbone", "heavy", "all"], default="heavy",
                        help="Atom selection mode for RMSF calculation (default: heavy atoms)")

    args = parser.parse_args()
    pdb_id = args.pdb_id
    mode = args.mode
    
    #if not pdb_id.isalnum() or len(pdb_id) != 4:
    #    print("Error: PDB ID is wrong >> " + pdb_id)
    #    sys.exit(1)

    get_rmsf_pdb(pdb_id, mode=mode)