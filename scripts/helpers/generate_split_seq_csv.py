import os
import csv
import argparse
from pathlib import Path

def generate_csv(split_name):
    """
    Generate a CSV file containing PDB IDs and their sequences from a split file.
    
    Args:
        split_name (str): Name of the split file (without .txt extension)
    """
    split_file = Path(f"./splits/{split_name}.txt")
    
    if not split_file.exists():
        print(f"Error: Split file {split_file} does not exist.")
        return
    
    output_csv = f"./splits/{split_name}_seqs.csv"
    
    with split_file.open('r') as f:
        pdb_ids = [line.strip() for line in f if line.strip()]
    
    with open(output_csv, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['name', 'seqres'])
        
        processed = 0
        for pdb_id in pdb_ids:
            seq_file = Path(f"./PDBs/{pdb_id}/{pdb_id}_seq.txt")
            
            if seq_file.exists():
                with seq_file.open('r') as f:
                    sequence = f.read().strip()
                
                csvwriter.writerow([pdb_id, sequence])
                processed += 1
            else:
                print(f"Warning: Sequence file not found for {pdb_id} at {seq_file}")
    
    print(f"CSV file generated: {output_csv}")
    print(f"Processed {processed} of {len(pdb_ids)} entries")

def main():
    parser = argparse.ArgumentParser(description='Generate CSV file with PDB IDs and sequences.')
    parser.add_argument('split_name', type=str, help='Name of the split file (without .txt extension)')
    
    args = parser.parse_args()
    generate_csv(args.split_name)

if __name__ == '__main__':
    main()