# Ensemble Benchmark
Comparing the accuracy of conformation predictors against experimental data with both global and local metrics.

## General Info

### Included Conformation Prediction Models
- Bioemu
- SAM2
- AlphaFlow
- Boltz2
- OpenFold (single conformation)

### CUDA Requirement
- The predictors are configured to use CUDA by default. If you prefer to run them on a CPU instead, you may need to manually update some dependencies and scripts (e.g., change `device='cuda:0'` to `device='cpu'`).
- The current dependencies are set up for an NVIDIA A800 GPU. If you're using a different GPU architecture, you may need to install a compatible version of PyTorch and CUDA.

## Installing dependencies

#### General Installation
- Clone this repository
- Clone the SAM2 repository: `./models/sam2`, and then install it as a package
- Clone the AlphaFlow repository: `./models/alphaflow`
- Clone and install OpenFold, and specify path in `./scripts/prepare_openfold.sh` and `./scripts/models/run_openfold.sh`. Alignment DBs are required, you can modify these paths in the scripts.
- Install the conda envs, etc. (Should install bioemu and boltz2)

#### Downloading Weights
- Download AlphaFlow weights with `./scripts/download_alphaflow_weights.sh`
- Other models will download weights on first run.


## Using Ensemble-Benchmark

### Preparing your Dataset
1. Create a [dataset_name].txt file in `./splits/`
2. In each line of the file should be a single four-character Protein Data Bank ID
3. Run `./scripts/prepare_split.sh ./splits/[dataset_name].txt` to prepare your dataset. This will download PDB, CIF, and MTZ files from PDBREDO to the PDB Storage Structure (see end).
4. **[AlphaFlow]** Run `./scripts/make_alphaflow_alignments.sh <dataset_name>` to make `./bin/alignment/*/a3m/*.a3m` alignment files for each PDB in the dataset. Runtime is around 20 seconds per PDB.
5. **[OpenFold]** Run `./scripts/prepare_openfold.sh <dataset_name>` to generate alignments for openfold. ⚠️ Runtime is around a few hours per PDB.
6. **[OpenFold]** Run `python ./scripts/helpers/download_alignment_cifs.py` to download CIF files for alignments prepared in the previous step. ⚠️ This will download 50 CIFs per PDB, and you may need to increase this number in the script depending on if OpenFold fails to find any CIF file. There are 500 total CIFs per Alignment, but for efficiency we only download 50.


### Inferencing Models
 The scripts below will deposit a multi-model PDB file (aligned with Phenix Phaser and mdtraj as a fallback) into the respective PDB folder, which will then be ready for analysis.

1. **Running BioEmu:**
   ```bash
   sh scripts/models/run_bioemu.sh <pdb_id> <samples>
   sh scripts/models/run_bioemu.sh 7lfo 1 # Example
   ```
   Output PDB: `./PDBs/*/*_bioemu.pdb`.
   <br>
   Misc Bin: `./PDBs/*/bioemu_bin/`

2. **Running SAM2:**
   ```bash
   sh scripts/models/run_sam2.sh <pdb_id> <samples>
   sh scripts/models/run_sam2.sh 7lfo 1 # Example
   ```
   Output PDB: `./PDBs/*/*_sam2.pdb`.
   <br>
   Misc Bin: `./PDBs/*/sam2_bin/`

3. **Running AlphaFlow:**
   ```bash
    sh scripts/models/run_alphaflow.sh <dataset_name> <samples>
    sh scripts/models/run_alphaflow.sh dataset 1 # Example
   ```
   Output PDB: `./PDBs/*/*_alphaflow.pdb`.
   <br>
   Misc Bin: `./PDBs/*/alphaflow_bin/`

4. **Running Boltz2**
   ```bash
    sh scripts/models/run_boltz2.sh <pdb_id> <samples>
    sh scripts/models/run_boltz2.sh 7lfo 1 # Example
   ```
   Output PDB: `./PDBs/*/*_boltz2.pdb`.
   <br>
   Misc Bin: `./PDBs/*/boltz2_bin/`


4. **Running OpenFold**
   ```bash
    sh scripts/models/run_openfold.sh <dataset_name>
    sh scripts/models/run_openfold.sh dataset # Example
    # Note: OpenFold will not terminate after finishing. You will need to Ctrl + C it when it finishes so the bash script can continue.
    # OpenFold only supports single conformation prediction. You will need to rename the Output PDB (to prevent overwriting) and run the above again. Repeat until you have the desired number of conformations. After that, use mdtraj to concatenate all the frames you just generated into one PDB and rename that to the final Output PDB file path for use in analysis. TODO: Turn this into a script
   ```
   Output PDB: `./PDBs/*/*_openfold.pdb`.
   <br>
   Misc Bin: `./PDBs/*/openfold_bin/`


### Obtaining Metrics
Now that all your PDBs are inferenced, analysis scripts can be ran to get metrics from them.

1. **R-free Calculation** - SFCalculator is used to determine the R-free value of each predicted PDB. Run it here:

    ```
      python ./scripts/analysis/get_rfrees.py <pdb_id> [--threads n]
    ```
    Output CSV Table `[predictor|rfree]`: `./PDBs/*/analysis/rfrees.csv`
    
2. **RMSF Calculation** - Root Mean Square Fluctuation Metric for each predicted PDB. This will be zero for OpenFold. Run it here:
    ```
      python ./scripts/analysis/get_rmsf.py <pdb_id>
    ```
    Output CSV Table `[predictor,residue,residue_aa,rmsf]`: `./PDBs/*/analysis/rmsf.csv`
    
3. **Density Fitness Metrics** - Local Metrics from density-fitness for each predicted PDB. Run it here:
    ```
      python ./scripts/analysis/get_density_fitness.py <pdb_id> [--threads n]
    ```
    Output CSV Table `[predictor,residue,JSON [ density_fitness_residue_obj ] ]`: `./PDBs/*/analysis/density_fitness.csv` where each conformation's residue has an object in the JSON array. The # of conformations generated will match the object count in the array.
    
4. **Secondary Structure** - Get the secondary structure of each residue. Run it here:
    ```
      python ./scripts/analysis/get_secondary_structure.py <pdb_id>
    ```
    Output CSV Table `[residue,secondary_structure]`: `./PDBs/*/analysis/secondary_structure.csv` where each secondary_structure is a one-character key for the secondary structure
    

### Visualizing and Summary
- various graphs in ./scripts/graphing . 
- summary stats and csv outputs in ./scripts/summary/

TODO: add docs for graphs and summary


--- 

### 📂 PDB Storage Structure

The `./PDBs/` directory is automatically created when running `prepare_split.sh`. Each PDB entry follows the general organization below. There may be additional files depending on the model and analysis performed.

```py
PDBs/
└── {pdb_id}/
    # Post prepare_split.sh
    ├── {pdb_id}_final.pdb       # PDBREDO Download
    ├── {pdb_id}_final.mtz       # PDBREDO Download
    ├── {pdb_id}_final.cif       # PDBREDO Download
    ├── {pdb_id}_seq.txt         # Amino Acid Sequence String
    ├── {pdb_id}_nowat.pdb       # PDB with only protein atoms

    # Post model inferencing
    ├── {pdb_id}_bioemu.pdb      # BIOEMU output with multi conformations
    ├── {pdb_id}_sam2.pdb        # SAM2 output with multi conformations
    ├── {pdb_id}_alphaflow.pdb   # AlphaFlow output with multi conformations
    ├── {pdb_id}_boltz2.pdb      # Boltz2 output with multi conformations
    ├── {pdb_id}_openfold.pdb    # OpenFold output with single conformation

    # Post analysis
    ├── analysis/
    │   ├── rfrees.csv           # R-free values for each predictor
    │   ├── rmsf.csv             # RMSF values for each predictor
    │   ├── density_fitness.csv   # Density fitness metrics for each predictor
    │   ├── secondary_structure.csv # Secondary structure for each residue
    │   └── ...                  # Any other analysis outputs

    # Model specific directories and dump
    ├── bioemu_bin/              # Misc files from BioEmu
    ├── sam2_bin/                # Misc files from SAM2
    ├── alphaflow_bin/           # Misc files from AlphaFlow
    ├── boltz2_bin/              # Misc files from Boltz2
    └── openfold_bin/            # Misc files from OpenFold

```

