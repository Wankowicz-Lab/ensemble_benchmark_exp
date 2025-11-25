# MMTBX Rfree calculation script

from __future__ import absolute_import, division, print_function

import argparse
import os
import sys
import csv

# CCTBX rfree script sent by Stephanie
# results: significantly decreased Rfree of pymol/mdtraj alignment methods

import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
from iotbx import reflection_file_reader
import libtbx.load_env  # makes sure cctbx env is set up

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def find_arrays(miller_arrays):
    """Pick out F_obs and R-free flags from the MTZ."""
    f_obs = None
    r_free_flags = None

    print("Available MTZ labels:")
    for ma in miller_arrays:
        label = ma.info().label_string()
        print("  ", label)

        # Try common Fobs labels
        if label in ("FOBS_X,SIGFOBS_X", "FP,SIGFP"):
            f_obs = ma

        # Try common R-free labels
        if label in ("R-free-flags", "FREE"):
            r_free_flags = ma

    return f_obs, r_free_flags


def get_rfree(pdb_file, mtz_file):
    pdb_file = pdb_file
    mtz_file = mtz_file
    free_flag_value = 0

    # Basic sanity checks
    if not os.path.isfile(pdb_file):
        sys.exit(f"ERROR: PDB file not found: {pdb_file}")
    if not os.path.isfile(mtz_file):
        sys.exit(f"ERROR: MTZ file not found: {mtz_file}")

    print(f"PDB file: {pdb_file}")
    print(f"MTZ file: {mtz_file}")
    print(f"Working directory: {os.getcwd()}")

    # Read PDB and build model
    pdb_inp = iotbx.pdb.input(file_name=pdb_file)
    model = mmtbx.model.manager(model_input=pdb_inp)

    # Read MTZ -> Miller arrays
    hkl_inp = reflection_file_reader.any_reflection_file(file_name=mtz_file)
    miller_arrays = hkl_inp.as_miller_arrays()

    f_obs, r_free_flags = find_arrays(miller_arrays)

    if f_obs is None:
        sys.exit(
            "ERROR: Could not find F_obs array in MTZ. "
            "Expected one of: 'FOBS_X,SIGFOBS_X' or 'FP,SIGFP'."
        )
    if r_free_flags is None:
        sys.exit(
            "ERROR: Could not find R-free flags array in MTZ. "
            "Expected one of: 'R-free-flags' or 'FREE'."
        )

    # Put on common set of Miller indices
    f_obs, r_free_flags = f_obs.common_sets(r_free_flags)

    # Convert to boolean mask: True = free-set reflections
    data = r_free_flags.data()
    free_mask = (data == free_flag_value)
    r_free_flags = r_free_flags.array(data=free_mask)

    print(
        "R-free flags counts -> free: %d, work: %d"
        % (r_free_flags.data().count(True), r_free_flags.data().count(False))
    )

    # Build f_model manager and compute R-values
    fmodel = mmtbx.f_model.manager(
        f_obs=f_obs,
        r_free_flags=r_free_flags,
        xray_structure=model.get_xray_structure(),
    )
    fmodel.update_all_scales()

    r_work = fmodel.r_work()
    r_free = fmodel.r_free()

    print("Final R-values:")
    print("  R_work = %6.4f" % r_work)
    print("  R_free = %6.4f" % r_free)

    return r_free



# original heatmap creation vvvvvvvvvvvvvvvv
# adjusted to use cctbx get_rfree method instead of SFCalculator 

def process_protein(protein_name, input_dir, dir_name):
    print(f"\n{'='*60}")
    print(f"Processing: {protein_name}")
    print(f"{'='*60}")
    
    mtz_file = None
    for mtz_path in Path(input_dir).glob(f"{protein_name}_final.mtz"):
        mtz_file = str(mtz_path)
        break
    
    if not mtz_file:
        print(f"ERROR: No MTZ file found for {protein_name}")
        return None
    
    print(f"Using MTZ: {mtz_file}")
    
    alignment_methods = {
        'MDTraj': 'mdtraj.pdb',
        'Phaser MR': 'phaser.pdb',
        'PyMOL CEAlign': 'pymol.pdb',
        'Refined': 'phenix_refined.pdb',
        'Refined (SA)': 'phenix_refined_annealing.pdb',
        'Deposited': 'deposited.pdb'
    }
    
    results = {'protein': dir_name}
    
    for method, filename in alignment_methods.items():
        pdb_path = os.path.join(input_dir, filename)
        
        if not os.path.exists(pdb_path):
            print(f"\n{method}: File not found - {pdb_path}")
            results[method] = None
            continue
        
        print(f"\n{method}: {pdb_path}")
        rfree = get_rfree(pdb_path, mtz_file)
        results[method] = rfree
    
    return results

def create_heatmap(df, output_path):
    heatmap_data = df.set_index('protein').T
    
    plt.figure(figsize=(max(10, len(df) * 1.5), 8))
    
    sns.heatmap(heatmap_data, 
                annot=True, 
                fmt='.4f', 
                cmap='RdYlGn_r',
                cbar_kws={'label': 'R-free'},
                vmin=0.0,
                vmax=0.6,
                linewidths=0.5,
                linecolor='gray')
    
    #plt.title('R-free Comparison Alignment Methods', fontsize=14, fontweight='bold')
    plt.xlabel('PDB ID', fontsize=24)
    plt.ylabel('Alignment Method', fontsize=24)
    plt.tight_layout()
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nHeatmap saved to: {output_path}")
    plt.close()

def main():
    if len(sys.argv) < 2:
        print("Usage: python save_rfrees.py <input_dir1> [input_dir2] ...")
        print("Example: python save_rfrees.py ./input/1aho.sam2 ./input/3k0n.openfold")
        sys.exit(1)
    
    input_dirs = sys.argv[1:]
    all_results = []
    
    for input_dir in input_dirs:
        if not os.path.exists(input_dir):
            print(f"WARNING: Directory not found: {input_dir}")
            continue
        
        dir_name = os.path.basename(input_dir.rstrip('/'))
        protein_name = dir_name.split('.')[0]
        
        results = process_protein(protein_name, input_dir, dir_name)
        if results:
            all_results.append(results)
            
            output_dir = f"./output/{dir_name}"
            os.makedirs(output_dir, exist_ok=True)
            
            df_single = pd.DataFrame([results])
            csv_path = os.path.join(output_dir, 'rfree_comparison_cctbx.csv')
            df_single.to_csv(csv_path, index=False)
            print(f"\nSaved results to: {csv_path}")
    
    if not all_results:
        print("\nNo results to process!")
        sys.exit(1)
    
    df_combined = pd.DataFrame(all_results)
    
    os.makedirs('./output', exist_ok=True)
    combined_csv = './output/all_rfree_comparison_cctbx.csv'
    df_combined.to_csv(combined_csv, index=False)
    
    heatmap_path = './output/rfree_heatmap_cctbx.png'
    create_heatmap(df_combined, heatmap_path)


if __name__ == '__main__':
    main()