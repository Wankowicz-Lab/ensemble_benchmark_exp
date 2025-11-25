# get_density_fitness.py <pdb_id> [--threads N]
# adds a density_fitness.csv to {PDB}/analysis as predictor|frame|metrics*
#   *metrics is a JSON in TXT format, not parsed.


import sys
import os
import json
import numpy as np
import mdtraj as md
import tempfile
import pandas as pd
import argparse
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor
import threading
import subprocess
import shutil

local_temp = threading.local()

# density-fitness rejecting PDBs for not having HEADER at top
def put_header(pdb_path):
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    if (lines[0].startswith("MODEL")):
        lines.remove(lines[0])  # Remove MODEL line if it exists, it causes density-fitness to error.

    if not lines or not lines[0].startswith("HEADER"):
        header = f"HEADER    {os.path.basename(pdb_path)}\n"
        lines.insert(0, header)

    for i in range(len(lines)):
        if lines[i].startswith("MODEL"):
            lines[i] = ""
    
    with open(pdb_path, 'w') as f:
        f.writelines(lines)


# spawns the internal density-fitness script and gets out a JSON/txt with local metrics
def get_density_fitness(pdb_path, mtz_path, out_path):
    local_metric_script = os.path.join(os.path.dirname(__file__), 'internal/get_local_metrics.sh')
    try: 
        result = subprocess.run(
            ["/bin/bash", local_metric_script, mtz_path, pdb_path, out_path],
            capture_output=True,
            text=True,
        )
    except Exception as e:
        print('[get_density_fitness.py] SPAWN DENSITY-FITNESS FAIL: ', e)
    
    if (result.stderr):
        print(f"[get_density_fitness.py] SPAWN DENSITY-FITNESS ERROR: {result.stderr.strip()}")
        return None

    if (os.path.exists(out_path)):
        with open(out_path, 'r') as f:
            txt = f.read()
            return txt.strip()
    
    return None



# This is the thread function for a given frame
def process_frame(frame, frame_idx, predictor, pdb_id, mtz_path, main_temp_dir):
    if not hasattr(local_temp, 'dir'):
        local_temp.dir = tempfile.mkdtemp(dir=main_temp_dir)
        
    frame_path = os.path.join(local_temp.dir, f"{pdb_id}_{predictor}_{frame_idx}.pdb")
    frame.save_pdb(frame_path)

    put_header(frame_path)

    try:
        metrics = get_density_fitness(frame_path, mtz_path, os.path.join(local_temp.dir, f"{pdb_id}_{predictor}_{frame_idx}_local_metrics.json"))
        
        if metrics is None:
            print(f"[get_density_fitness.py] No metrics returned for {predictor} frame {frame_idx}")
            return None
        
        print(f"[get_density_fitness.py] CALCULATED [{predictor} - {frame_idx}] Local Metrics")
        return {
            "predictor": predictor,
            "frame": frame_idx,
            "metrics": json.loads(metrics),
        }
    except Exception as e:
        print(f"[get_density_fitness.py] ERROR processing {predictor} frame {frame_idx}: {e}")
        return None


def make_density_fitnesses(pdb_id, num_threads=5):
    predictors = ['bioemu', 'alphaflow', 'sam2', 'boltz2', 'openfold']
    MTZ_PATH = f"./PDBs/{pdb_id}/{pdb_id.lower()}_final.mtz"
    metrics = []
    
    with tempfile.TemporaryDirectory() as main_temp_dir:
        for i, predictor in enumerate(predictors):
            ensemble_path = f"./PDBs/{pdb_id}/{pdb_id.lower()}_{predictor}.pdb"
            
            if not os.path.exists(ensemble_path):
                print(f"[get_density_fitness.py] Warning: {ensemble_path} not found, skipping...")
                continue
                
            try:
                frames = md.load(ensemble_path)
                print(f"[get_density_fitness.py] Processing {predictor} with {frames.n_frames} frames using {num_threads} threads")
                
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
                        if result is not None:
                            metrics.append(result)
                        
                        completed += 1
                        if completed % num_threads == 0 or completed == len(futures):
                            print(f"[get_density_fitness.py] Progress: {completed}/{len(futures)} frames processed for {predictor}")
                            
                print(f"[get_density_fitness.py] Completed {predictor} for {pdb_id}! {len(predictors)-i-1} predictors left!")
                
            except Exception as e:
                print(f"[get_density_fitness.py] ERROR processing {predictor}: {e}")
                continue

    csv = pd.DataFrame([r for r in metrics if r["metrics"] is not None])
    if not csv.empty:
        os.makedirs(f"./PDBs/{pdb_id}/analysis", exist_ok=True)
        csv.to_json(f"./PDBs/{pdb_id}/analysis/density_fitness.json", index=False)
        print(f"[get_density_fitness.py] Saved density_fitness.json for {pdb_id} in ./PDBs/{pdb_id}/analysis/density_fitness.csv")
    else:
        print(f"[get_density_fitness.py] No valid Density Fitness values calculated")




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Density Fitness values for protein model ensembles by predictors")
    parser.add_argument("pdb_id", help="PDB ID to process")
    parser.add_argument("--threads", "-t", type=int, default=5, help="Number of threads to use")
    
    args = parser.parse_args()
    pdb_id = args.pdb_id
    
    if not pdb_id.isalnum() or len(pdb_id) != 4:
        print("Error: PDB ID is wrong >> " + pdb_id)
        sys.exit(1)

    make_density_fitnesses(pdb_id, args.threads)