import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


from SFC_Torch import SFcalculator

def get_rfree(pdb_file, mtz_file):
    sfcalculator = SFcalculator(pdb_file, mtz_file, expcolumns=['FP', 'SIGFP'], set_experiment=True, freeflag='FREE', testset_value=0)
    sfcalculator.inspect_data(verbose=False) 
    sfcalculator.calc_fprotein(atoms_position_tensor=None, atoms_biso_tensor=None, atoms_occ_tensor=None, atoms_aniso_uw_tensor=None)
    sfcalculator.calc_fsolvent()
    sfcalculator.init_scales(requires_grad=True)
    Fmodel = sfcalculator.calc_ftotal()
    return sfcalculator.r_free.item()



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
    
    #plt.title('R-free Comparison: Alignment Methods', fontsize=14, fontweight='bold')
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
            csv_path = os.path.join(output_dir, 'rfree_comparison.csv')
            df_single.to_csv(csv_path, index=False)
            print(f"\nSaved results to: {csv_path}")
    
    if not all_results:
        print("\nNo results to process!")
        sys.exit(1)
    
    df_combined = pd.DataFrame(all_results)
    
    os.makedirs('./output', exist_ok=True)
    combined_csv = './output/all_rfree_comparison.csv'
    df_combined.to_csv(combined_csv, index=False)
    
    heatmap_path = './output/rfree_heatmap.png'
    create_heatmap(df_combined, heatmap_path)


if __name__ == '__main__':
    main()