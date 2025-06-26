# python pdb.py <pdb> <predictor>


import pandas as pd
import numpy as np
import os
import sys
import json
import argparse

def load_rmsf_data(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/rmsf.csv"
    
    if not os.path.exists(file_path):
        print(f"[pdb.py] RMSF file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        print(f"[pdb.py] Error reading RMSF data: {e}")
        return pd.DataFrame()


def load_ediam_data(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/density_fitness.csv"
    
    if not os.path.exists(file_path):
        print(f"[pdb.py] Density fitness file not found: {file_path}")
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
                    print(f"[pdb.py] Error parsing metrics: {e}")
            
            ediam_df = pd.DataFrame(all_data)
            
            # More Statistics per stat:
            
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
        print(f"[pdb.py] Error loading density fitness data: {e}")
        return pd.DataFrame()


def load_secondary_structure(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/secondary_structure.csv"
    
    if not os.path.exists(file_path):
        print(f"[pdb.py] Secondary structure file not found: {file_path}")
        return pd.DataFrame()
    
    try:
        return pd.read_csv(file_path)
    except Exception as e:
        print(f"[pdb.py] Error reading secondary structure data: {e}")
        return pd.DataFrame()

def load_rfree(pdb_id):
    file_path = f"./PDBs/{pdb_id}/analysis/rfrees.csv"
    
    if not os.path.exists(file_path):
        print(f"[pdb.py] R-free file not found: {file_path}")
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

if __name__ == "__main__":
    print("Usage: python pdb.py <pdb_id> [--predictor <predictor_name>] [--output_csv <output_file>]")

    argparse = argparse.ArgumentParser(description="Output information about predictions from PDB ID.")
    argparse.add_argument('pdb_id', type=str, help='PDB ID to process')
    argparse.add_argument('--predictor', type=str, help='Predictor name to filter data')
    argparse.add_argument('--output_csv', type=str, help='Output CSV file to save results')

    args = argparse.parse_args()
    pdb_id = args.pdb_id


    if not pdb_id:
        print("Invalid Usage")
        sys.exit(1)


    rmsf_data = load_rmsf_data(pdb_id)
    ediam_data = load_ediam_data(pdb_id)
    ss_data = load_secondary_structure(pdb_id)
    ss_data['residue'] = ss_data['residue'].astype(int) - 1 # residues are 0-based here


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

    predictorFilter = args.predictor
    if predictorFilter:
        merged_data = merged_data[merged_data['predictor'] == predictorFilter]

    if args.output_csv:
        output_file = args.output_csv
        merged_data.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    
    merged_data = merged_data.sort_values(by=['residue', 'predictor']).reset_index(drop=True)


    residue_count = merged_data['residue'].nunique()

    r_frees = load_rfree(pdb_id)
    rfree_sam2 = r_frees['sam2'] if r_frees else None
    rfree_alphaflow = r_frees['alphaflow'] if r_frees else None
    rfree_bioemu = r_frees['bioemu'] if r_frees else None

    print("\nPDB Statistics (PDB ID:", pdb_id.upper() + ")")

    print("\n====== PREDICTED ENSEMBLE R-FREES ======")
    print("SAM2: ", rfree_sam2)
    print("ALPHAFLOW: ", rfree_alphaflow)
    print("BIOEMU:", rfree_bioemu)

    for _, thisPredictor in enumerate(['sam2', 'alphaflow', 'bioemu']):
        if (args.predictor and args.predictor != thisPredictor):
            continue
        print(f"\n====== {thisPredictor.upper()} - DENSITY FITNESS STATISTICS ======")
        thisMergedData = merged_data[merged_data['predictor'] == thisPredictor]
        avg_esdiam = thisMergedData['EDIAm_mean'].mean()
        sd_esdiam = thisMergedData['EDIAm_std'].mean()
        iqr_esdiam = thisMergedData['EDIAm_iqr'].mean()
        min_esdiam = thisMergedData['EDIAm_min'].min()
        max_esdiam = thisMergedData['EDIAm_max'].max()

        avg_rsccs = thisMergedData['RSCCS_mean'].mean()
        sd_rsccs = thisMergedData['RSCCS_std'].mean()
        iqr_rsccs = thisMergedData['RSCCS_iqr'].mean()
        min_rsccs = thisMergedData['RSCCS_min'].min()
        max_rsccs = thisMergedData['RSCCS_max'].max()
        
        avg_rsr = thisMergedData['RSR_mean'].mean()
        sd_rsr = thisMergedData['RSR_std'].mean()
        iqr_rsr = thisMergedData['RSR_iqr'].mean()
        min_rsr = thisMergedData['RSR_min'].min()
        max_rsr = thisMergedData['RSR_max'].max()

        avg_rmsf = thisMergedData['RMSF'].mean()
        sd_rmsf = thisMergedData['RMSF'].std()
        iqr_rmsf = thisMergedData['RMSF'] .quantile(0.75) - thisMergedData['RMSF'].quantile(0.25)
        min_rmsf = thisMergedData['RMSF'].min()
        max_rmsf = thisMergedData['RMSF'].max()

            
        avg_stat_df = pd.DataFrame(
            columns=['stat', 'mean', 'std', 'iqr', 'min', 'max'],
            data=[
                ['RMSF', avg_rmsf, sd_rmsf, iqr_rmsf, min_rmsf, max_rmsf],
                ['EDIAm', avg_esdiam, sd_esdiam, iqr_esdiam, min_esdiam, max_esdiam],
                ['RSCCS', avg_rsccs, sd_rsccs, iqr_rsccs, min_rsccs, max_rsccs],
                ['RSR', avg_rsr, sd_rsr, iqr_rsr, min_rsr, max_rsr]
            ]

        )
        print(avg_stat_df)

    print("\n====== PER-RESIDUE STATS ======")
    print("Residue Count:", residue_count)

    print(merged_data)


