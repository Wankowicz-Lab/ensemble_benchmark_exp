import mdtraj as md
import os
import sys

def align_with_mdtraj(ensemble_path, topology_path):
    # Load the reference structure (topology)
    topology = md.load(topology_path)
    # Load the ensemble (multi-frame PDB or trajectory)
    ensemble = md.load(ensemble_path)
    # Choose atom indices for alignment (must be common to both)
    atom_indices = topology.topology.select("name CA")  # You can change to "backbone" or custom selection
    # Align all frames in the ensemble to the topology
    ensemble.superpose(topology, atom_indices=atom_indices)
    # Save the aligned ensemble if needed
    ensemble.save(ensemble_path)
    print(f'[align_with_mdtraj.py] Aligned ensemble to topology and saved to {ensemble_path}')

if __name__ == "__main__":
    ensemble_path = sys.argv[1]
    topology_path = sys.argv[2]

    if not os.path.exists(ensemble_path):
        print(f"Error: DCD file '{ensemble_path}' does not exist.")
        sys.exit(1)
    if not os.path.exists(topology_path):
        print(f"Error: Topology file '{topology_path}' does not exist.")
        sys.exit(1)

    
    align_with_mdtraj(ensemble_path, topology_path)
    sys.exit(0)