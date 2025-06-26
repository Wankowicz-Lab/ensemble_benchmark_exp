# python dataset.py <dataset> [--predictor <predictor>] [--output_csv <output_file>]


import pandas as pd
import numpy as np
import os
import sys
import json
import argparse

def load_rmsf_data(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/rmsf.csv"
    
    if not os.path.exists(file_path):
        print(f"[dataset.py] RMSF file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        print(f"[dataset.py] Error reading RMSF data: {e}")
        return pd.DataFrame()


def load_ediam_data(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/density_fitness.csv"
    
    if not os.path.exists(file_path):
        print(f"[dataset.py] Density fitness file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        lines = []
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    lines.append({
                        'predictor': line.split(',')[0].strip(),
                        'frame': line.split(',')[1].strip(),
                        'metrics': (','.join(line.split(',')[2:]).strip())[1:-1]  # Remove surrounding quotes cuz JSON parsing
                    })
        
        df = pd.DataFrame(lines)
        all_data = []
        
        if 'metrics' in df.columns:
            for _, row in df.iterrows():
                predictor = row['predictor']
                frame = row['frame']
                if (_ == 0):
                    continue 

                try:
                    metrics = json.loads(row['metrics'])
                    
                    for residue_data in metrics:
                        if 'EDIAm' in residue_data and 'seqID' in residue_data:
                            all_data.append({
                                'predictor': predictor,
                                'frame': frame,
                                'residue': (residue_data['pdb']['seqNum'] - 1) if predictor == 'alphaflow' else residue_data['pdb']['seqNum'], # alphaflow starts residue number at 1 
                                'EDIAm': residue_data['EDIAm'],
                                'RSCCS': residue_data['RSCCS'],
                                'RSR': residue_data['RSR'],
                                'SRSR': residue_data['SRSR'],
                                'aa': residue_data.get('compID', '')
                            })
                except Exception as e:
                    print(f"[dataset.py] Error parsing metrics: {e}")
            
            ediam_df = pd.DataFrame(all_data)
            
            
            stats_ediam = ediam_df.groupby(['predictor', 'residue', 'aa'])['EDIAm'].agg([
                'mean', 'std', 'min', 'max',
                lambda x: x.quantile(0.25),  # Q1
                lambda x: x.quantile(0.75)   # Q3
            ]).reset_index()
            stats_ediam.columns = ['predictor', 'residue', 'aa', 'EDIAm_mean', 'EDIAm_std', 'EDIAm_min', 'EDIAm_max', 'EDIAm_q1', 'EDIAm_q3']
            stats_ediam['EDIAm_iqr'] = stats_ediam['EDIAm_q3'] - stats_ediam['EDIAm_q1']
            
            stats_rsccs = ediam_df.groupby(['predictor', 'residue', 'aa'])['RSCCS'].agg([
                'mean', 'std', 'min', 'max',
                lambda x: x.quantile(0.25),  # Q1
                lambda x: x.quantile(0.75)   # Q3
            ]).reset_index()
            stats_rsccs.columns = ['predictor', 'residue', 'aa', 'RSCCS_mean', 'RSCCS_std', 'RSCCS_min', 'RSCCS_max', 'RSCCS_q1', 'RSCCS_q3']
            stats_rsccs['RSCCS_iqr'] = stats_rsccs['RSCCS_q3'] - stats_rsccs['RSCCS_q1']
            
            stats_rsr = ediam_df.groupby(['predictor', 'residue', 'aa'])['RSR'].agg([
                'mean', 'std', 'min', 'max',
                lambda x: x.quantile(0.25),  # Q1
                lambda x: x.quantile(0.75)   # Q3
            ]).reset_index()
            stats_rsr.columns = ['predictor', 'residue', 'aa', 'RSR_mean', 'RSR_std', 'RSR_min', 'RSR_max', 'RSR_q1', 'RSR_q3']
            stats_rsr['RSR_iqr'] = stats_rsr['RSR_q3'] - stats_rsr['RSR_q1']
            
            ediam_df = pd.merge(stats_ediam, stats_rsccs, on=['predictor', 'residue', 'aa'], how='outer')
            ediam_df = pd.merge(ediam_df, stats_rsr, on=['predictor', 'residue', 'aa'], how='outer')

            ediam_df = ediam_df[[
                'predictor', 'residue', 'aa', 
                'EDIAm_mean', 'EDIAm_std', 'EDIAm_iqr', 'EDIAm_min', 'EDIAm_max',
                'RSCCS_mean', 'RSCCS_std', 'RSCCS_iqr', 'RSCCS_min', 'RSCCS_max',
                'RSR_mean', 'RSR_std', 'RSR_iqr', 'RSR_min', 'RSR_max'
            ]]
            return ediam_df
        
        return pd.DataFrame()
        
    except Exception as e:
        print(f"[dataset.py] Error loading density fitness data: {e}")
        return pd.DataFrame()


def load_secondary_structure(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/secondary_structure.csv"
    
    if not os.path.exists(file_path):
        print(f"[dataset.py] Secondary structure file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        print(f"[dataset.py] Error reading secondary structure data: {e}")
        return pd.DataFrame()

def load_rfree(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/rfrees.csv"
    
    if not os.path.exists(file_path):
        print(f"[dataset.py] R-free file not found: {file_path}")
        return None
    
    rfree_df = pd.read_csv(file_path)
    
    sam2_df = rfree_df[rfree_df['predictor'] == 'sam2'].reset_index(drop=True)
    alphaflow_df = rfree_df[rfree_df['predictor'] == 'alphaflow'].reset_index(drop=True)
    bioemu_df = rfree_df[rfree_df['predictor'] == 'bioemu'].reset_index(drop=True)

    return {
        'sam2': sam2_df['rfree'].get(0, None),
        'alphaflow': alphaflow_df['rfree'].get(0, None),
        'bioemu': bioemu_df['rfree'].get(0, None)
    }

def load_pdb_list(split_name):
    split_file = f"./splits/{split_name}.txt"
    
    if not os.path.exists(split_file):
        print(f"[dataset.py] Split file not found: {split_file}")
        return []
    
    try:
        with open(split_file, 'r') as f:
            pdb_ids = [line.strip() for line in f if line.strip()]
        return pdb_ids
    except Exception as e:
        print(f"[dataset.py] Error reading split file: {e}")
        return []

def process_pdb_stats(pdb_id):
    rmsf_data = load_rmsf_data(pdb_id)
    ediam_data = load_ediam_data(pdb_id)
    ss_data = load_secondary_structure(pdb_id)
    r_frees = load_rfree(pdb_id)
    
    if ss_data.empty or rmsf_data.empty or ediam_data.empty:
        return None
    
    ss_data['residue'] = ss_data['residue'].astype(int) - 1  # residues are 0-based here

    merged_data = pd.merge(
        rmsf_data,
        ediam_data,
        on=['predictor', 'residue'],
        how='inner',
        suffixes=('_rmsf', '_ediam')
    )

    merged_data = pd.merge(
        merged_data,
        ss_data,
        on=['residue'],
        how='left'
    )

    merged_data['RMSF'] = merged_data['rmsf'].astype(float)
    merged_data = merged_data[[
        'residue', 'aa', 'secondary_structure', 'predictor', 'RMSF', 
        'EDIAm_mean', 'EDIAm_std', 'EDIAm_iqr', 'EDIAm_min', 'EDIAm_max',
        'RSCCS_mean', 'RSCCS_std', 'RSCCS_iqr', 'RSCCS_min', 'RSCCS_max',
        'RSR_mean', 'RSR_std', 'RSR_iqr', 'RSR_min', 'RSR_max'
    ]]

    pdb_stats = []
    
    for predictor in ['sam2', 'alphaflow', 'bioemu']:
        predictor_data = merged_data[merged_data['predictor'] == predictor]
        
        if predictor_data.empty:
            continue
            
        residue_count = predictor_data['residue'].nunique()
        
        stats = {
            'pdb_id': pdb_id,
            'predictor': predictor,
            'residue_count': residue_count,
            'rfree': r_frees[predictor] if r_frees else None,
            
            # RMSF stats
            'rmsf_mean': predictor_data['RMSF'].mean(),
            'rmsf_std': predictor_data['RMSF'].std(),
            'rmsf_min': predictor_data['RMSF'].min(),
            'rmsf_max': predictor_data['RMSF'].max(),
            'rmsf_iqr': predictor_data['RMSF'].quantile(0.75) - predictor_data['RMSF'].quantile(0.25),
            
            # EDIAm stats
            'ediam_mean': predictor_data['EDIAm_mean'].mean(),
            'ediam_std': predictor_data['EDIAm_std'].mean(),
            'ediam_min': predictor_data['EDIAm_min'].min(),
            'ediam_max': predictor_data['EDIAm_max'].max(),
            'ediam_iqr': predictor_data['EDIAm_iqr'].mean(),
            
            # RSCCS stats
            'rsccs_mean': predictor_data['RSCCS_mean'].mean(),
            'rsccs_std': predictor_data['RSCCS_std'].mean(),
            'rsccs_min': predictor_data['RSCCS_min'].min(),
            'rsccs_max': predictor_data['RSCCS_max'].max(),
            'rsccs_iqr': predictor_data['RSCCS_iqr'].mean(),
            
            # RSR stats
            'rsr_mean': predictor_data['RSR_mean'].mean(),
            'rsr_std': predictor_data['RSR_std'].mean(),
            'rsr_min': predictor_data['RSR_min'].min(),
            'rsr_max': predictor_data['RSR_max'].max(),
            'rsr_iqr': predictor_data['RSR_iqr'].mean(),
        }
        
        pdb_stats.append(stats)
    
    return pdb_stats

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate per-PDB statistics for a dataset split.")
    parser.add_argument('split_name', type=str, help='Name of the split (e.g., train, test, validation)')
    parser.add_argument('--predictor', type=str, help='Predictor name to filter data')
    parser.add_argument('--output_csv', type=str, help='Output CSV file to save results')

    args = parser.parse_args()
    split_name = args.split_name

    pdb_ids = load_pdb_list(split_name)
    
    if not pdb_ids:
        print(f"No PDB IDs found in split: {split_name}")
        sys.exit(1)

    print(f"Processing {len(pdb_ids)} PDBs from split: {split_name}")

    all_stats = []
    
    for i, pdb_id in enumerate(pdb_ids):
        print(f"Processing PDB {i+1}/{len(pdb_ids)}: {pdb_id}")
        
        pdb_stats = process_pdb_stats(pdb_id)
        if pdb_stats:
            all_stats.extend(pdb_stats)
        else:
            print(f"[dataset.py] Skipping {pdb_id} due to missing data")

    if not all_stats:
        print("No valid statistics found for any PDB in the split")
        sys.exit(1)

    stats_df = pd.DataFrame(all_stats)

    if args.predictor:
        stats_df = stats_df[stats_df['predictor'] == args.predictor]

    if args.output_csv:
        stats_df.to_csv(args.output_csv, index=False)
        print(f"Results saved to {args.output_csv}")

    # show summary statistics
    print(f"\n====== DATASET STATISTICS ({split_name.upper()}) ======")
    print(f"Total PDBs processed: {stats_df['pdb_id'].nunique()}")
    print(f"Total records: {len(stats_df)}")

    for predictor in ['sam2', 'alphaflow', 'bioemu']:
        if args.predictor and args.predictor != predictor:
            continue
            
        predictor_stats = stats_df[stats_df['predictor'] == predictor]
        
        if predictor_stats.empty:
            continue
            
        print(f"\n====== {predictor.upper()} SUMMARY ======")
        print(f"PDBs with {predictor} data: {len(predictor_stats)}")
        
        rfree_data = predictor_stats['rfree'].dropna()
      
        new_stats_df = pd.DataFrame(
            columns=['stat', 'mean', 'std', 'min', 'max'],
            data=[
                ['R-free', rfree_data.mean(), rfree_data.std(), rfree_data.min(), rfree_data.max()],
                ['RMSF', predictor_stats['rmsf_mean'].mean(), predictor_stats['rmsf_std'].mean(), predictor_stats['rmsf_min'].min(), predictor_stats['rmsf_max'].max()],
                ['EDIAm', predictor_stats['ediam_mean'].mean(), predictor_stats['ediam_std'].mean(), predictor_stats['ediam_min'].min(), predictor_stats['ediam_max'].max()],
                ['RSCCS', predictor_stats['rsccs_mean'].mean(), predictor_stats['rsccs_std'].mean(), predictor_stats['rsccs_min'].min(), predictor_stats['rsccs_max'].max()],
                ['RSR', predictor_stats['rsr_mean'].mean(), predictor_stats['rsr_std'].mean(), predictor_stats['rsr_min'].min(), predictor_stats['rsr_max'].max()]
            ]
        )
        print(new_stats_df)
        
    print(f"\n====== PER-PDB STATISTICS ======")
    print(stats_df.head())