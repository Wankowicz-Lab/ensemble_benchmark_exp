import mdtraj as md

ensemble_path = '/dors/wankowicz_lab/castelt/EnsembleBenchmark/PDBs/7lfo/7lfo_bioemu.pdb'
topology_path = '/dors/wankowicz_lab/castelt/EnsembleBenchmark/PDBs/7lfo/7lfo_final.pdb'
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
