import os
import sys
import mdtraj as md
 
def dcd_to_pdb(dcd_file, top_file, output_file):
    traj = md.load(dcd_file, top=top_file)
    traj.save_pdb(output_file)
    print(f'[dcd_to_pdb.py] Saved trajectory to {output_file}')

if __name__ == "__main__":
    dcd_path = sys.argv[1]
    top_path = sys.argv[2]
    output_path = sys.argv[3]

    if not os.path.exists(dcd_path):
        print(f"Error: DCD file '{dcd_path}' does not exist.")
        sys.exit(1)
    if not os.path.exists(top_path):
        print(f"Error: Topology file '{top_path}' does not exist.")
        sys.exit(1)
    if not output_path:
        print("Error: Output file path must be specified.")
        sys.exit(1)
    
    dcd_to_pdb(dcd_path, top_path, output_path)
    sys.exit(0)