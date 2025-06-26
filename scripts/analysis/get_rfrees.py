# get_rfrees.py <pdb_id> [--threads N]
# adds a rfrees.csv to {PDB}/analysis as predictor|rfree

import sys
import os
import json
import numpy as np
import mdtraj as md
import tempfile
from SFC_Torch import SFcalculator
import pandas as pd
import argparse
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import threading

local_temp = threading.local()

def get_rfree(pdb_file, mtz_file, predictor):
    sfcalculator = SFcalculator(pdb_file, mtz_file, expcolumns=['FP', 'SIGFP'], set_experiment=True, freeflag='FREE', testset_value=0)
    sfcalculator.inspect_data(verbose=False) 
    sfcalculator.calc_fprotein(atoms_position_tensor=None, atoms_biso_tensor=None, atoms_occ_tensor=None, atoms_aniso_uw_tensor=None)
    sfcalculator.calc_fsolvent()
    sfcalculator.init_scales(requires_grad=True)
    Fmodel = sfcalculator.calc_ftotal()
    return { "rfree": sfcalculator.r_free, "predictor": predictor }


def make_rfrees(pdb_id, num_threads=3):
    predictors = ['bioemu', 'alphaflow', 'sam2', 'boltz2', 'openfold']
    MTZ_PATH = f"./PDBs/{pdb_id}/{pdb_id.lower()}_final.mtz"
    rFrees = []
    
    with tempfile.TemporaryDirectory() as main_temp_dir:
        with ThreadPoolExecutor(max_workers=num_threads) as executor:
            futures = []
            for predictor in predictors:
                ensemble_path = f"./PDBs/{pdb_id}/{pdb_id.lower()}_{predictor}.pdb"
                
                if not os.path.exists(ensemble_path):
                    print(f"[get_rfrees.py] Warning: {ensemble_path} not found, skipping...")
                    continue
                
                future = executor.submit(get_rfree, ensemble_path, MTZ_PATH, predictor)
                futures.append(future)

            completed = 0
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                if result["rfree"] is not None:
                    rFrees.append({
                        "predictor": result["predictor"],
                        "rfree": float(result["rfree"]),
                    })
                        
                completed += 1
                print(f"[get_rfrees.py] Progress: {completed}/{len(futures)} predictors processed")
        

    csv = pd.DataFrame([r for r in rFrees if r["rfree"] is not None])
    if not csv.empty:
        os.makedirs(f"./PDBs/{pdb_id}/analysis", exist_ok=True)
        csv.to_csv(f"./PDBs/{pdb_id}/analysis/rfrees.csv", index=False)
        print(f"[get_rfrees.py] Saved rfrees.csv for {pdb_id} in ./PDBs/{pdb_id}/analysis/rfrees.csv")
    else:
        print(f"[get_rfrees.py] No valid R-free values calculated")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate R-free values for protein model ensembles by predictors")
    parser.add_argument("pdb_id", help="PDB ID to process")
    parser.add_argument("--threads", "-t", type=int, default=3, help="Number of threads to use")
    
    args = parser.parse_args()
    pdb_id = args.pdb_id
    
    if not pdb_id.isalnum() or len(pdb_id) != 4:
        print("Error: PDB ID is wrong >> " + pdb_id)
        sys.exit(1)

    make_rfrees(pdb_id, args.threads)

    