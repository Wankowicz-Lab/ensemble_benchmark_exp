# get_rfrees.py <pdb_id> [--threads N]
# adds a rfrees.csv to {PDB}/analysis as predictor|frame|rfree

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

def get_rfree(pdb_file, mtz_file):
    sfcalculator = SFcalculator(pdb_file, mtz_file, expcolumns=['FP', 'SIGFP'], set_experiment=True, freeflag='FREE', testset_value=0)
    sfcalculator.inspect_data(verbose=False) 
    sfcalculator.calc_fprotein(atoms_position_tensor=None, atoms_biso_tensor=None, atoms_occ_tensor=None, atoms_aniso_uw_tensor=None)
    sfcalculator.calc_fsolvent()
    sfcalculator.init_scales(requires_grad=True)
    Fmodel = sfcalculator.calc_ftotal()
    return sfcalculator.r_free

# This is the thread function for agiven frame
def process_frame(frame, frame_idx, predictor, pdb_id, mtz_path, main_temp_dir):
    if not hasattr(local_temp, 'dir'):
        local_temp.dir = tempfile.mkdtemp(dir=main_temp_dir)
        
    frame_path = os.path.join(local_temp.dir, f"{pdb_id}_{predictor}_{frame_idx}.pdb")
    frame.save_pdb(frame_path)
    
    try:
        rfree = get_rfree(frame_path, mtz_path)
        print(f"[get_rfrees.py] CALCULATED [{predictor} - {frame_idx}] Rfree: {rfree}")
        return {
            "predictor": predictor,
            "frame": frame_idx,
            "rfree": float(rfree)
        }
    except Exception as e:
        print(f"[get_rfrees.py] ERROR processing {predictor} frame {frame_idx}: {e}")
        return {
            "predictor": predictor,
            "frame": frame_idx,
            "rfree": None,
            "error": str(e)
        }

def make_rfrees(pdb_id, num_threads=4):
    predictors = ['bioemu', 'alphaflow', 'sam2', 'boltz2', 'openfold']
    MTZ_PATH = f"./PDBs/{pdb_id}/{pdb_id.lower()}_final.mtz"
    rFrees = []
    
    with tempfile.TemporaryDirectory() as main_temp_dir:
        for i, predictor in enumerate(predictors):
            ensemble_path = f"./PDBs/{pdb_id}/{pdb_id.lower()}_{predictor}.pdb"
            
            if not os.path.exists(ensemble_path):
                print(f"[get_rfrees.py] Warning: {ensemble_path} not found, skipping...")
                continue
                
            try:
                frames = md.load(ensemble_path)
                print(f"[get_rfrees.py] Processing {predictor} with {frames.n_frames} frames using {num_threads} threads")
                
                # Process frames in parallel
                with ThreadPoolExecutor(max_workers=num_threads) as executor:
                    futures = []
                    for frame_idx, frame in enumerate(frames):
                        future = executor.submit(
                            process_frame, 
                            frame, 
                            frame_idx, 
                            predictor, 
                            pdb_id, 
                            MTZ_PATH, 
                            main_temp_dir
                        )
                        futures.append(future)
                    
                    completed = 0
                    for future in concurrent.futures.as_completed(futures):
                        result = future.result()
                        if result["rfree"] is not None:
                            rFrees.append(result)
                        
                        completed += 1
                        if completed % 5 == 0 or completed == len(futures):
                            print(f"[get_rfrees.py] Progress: {completed}/{len(futures)} frames processed for {predictor}")
                            
                print(f"[get_rfrees.py] Completed {predictor} for {pdb_id}! {len(predictors)-i-1} predictors left!")
                
            except Exception as e:
                print(f"[get_rfrees.py] ERROR processing {predictor}: {e}")
                continue

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
    parser.add_argument("--threads", "-t", type=int, default=4, help="Number of threads to use")
    
    args = parser.parse_args()
    pdb_id = args.pdb_id
    
    if not pdb_id.isalnum() or len(pdb_id) != 4:
        print("Error: PDB ID is wrong >> " + pdb_id)
        sys.exit(1)

    make_rfrees(pdb_id, args.threads)

    