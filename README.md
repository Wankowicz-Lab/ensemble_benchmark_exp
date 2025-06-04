# Ensemble Benchmark
 Ensemble Benchmark aims to use **bioemu, sam2, and alphaflow** to predict multiple conformations given a PDB ID. The instructions below show how to get started.




## Running Instructions 

...for now 

### Installing dependencies
- Install the conda env, etc.
- Clone the SAM2 repository: 
  ```bash
  git clone <SAM_REPO_URL> ./models/sam2
  ```

---

## Using Ensemble-Benchmark

1. **Create a new `.txt` file in `./splits`**
   - The `.txt` file should contain a list of PDB IDs.
   
2. **Run the preparation script**
    
   This script will download the PDB files to ./PDBs and prepare the inputs needed for prediction models. 
   ```bash
   python scripts/prepare_split.py splits/dataset.txt
   ```
   (Replace `dataset.txt` with your split name)
   
   

3. **Run BioEmu Individually:**
   ```bash
   sh scripts/models/run_bioemu.sh <pdb_id> <samples>
   sh scripts/models/run_bioemu.sh 7lfo 1
   ```

4. **Run SAM2 Individually:**
   ```bash
   sh scripts/models/run_sam2.sh <pdb_id> <samples>
   sh scripts/models/run_sam2.sh 7lfo 1
   ```

**TODO:** The main handler should call all these scripts and then pipeline the analysis. This supports individual calling given a PDB at the moment.

---

## 📚 Additional Documentation

### 📂 PDB Storage Structure

The `./PDBs/` directory is automatically created when running `prepare_split.py`. Each PDB entry follows this organization:

```
PDBs/
└── {pdb_id}/
    ├── {pdb_id}_final.pdb       # PDBREDO or Original Download
    ├── {pdb_id}_final.mtz       # Only present when PDB available from PDBREDO
    ├── {pdb_id}_sequence.txt    # Amino Acid Sequence String
    ├── {pdb_id}_nowat.pdb       # PDB with water and MG ions removed
    ├── {pdb_id}_bioemu.pdb      # BIOEMU output with multi conformations
    ├── {pdb_id}_sam2.pdb        # SAM2 output with multi conformations
    ├── sam2.top.pdb             # SAM2 Topology PDB File
    └── sam2.traj.dcd            # SAM2 Trajectory DCD File
```

#### File Details:

<table>
  <tr>
    <th colspan="2">📥 Initial Files</th>
  </tr>
  <tr>
    <td><code>{pdb_id}_final.pdb</code></td>
    <td>PDBREDO or Original Download structure</td>
  </tr>
  <tr>
    <td><code>{pdb_id}_final.mtz</code></td>
    <td>Only present when PDB available from PDBREDO</td>
  </tr>
  <tr>
    <td><code>{pdb_id}_sequence.txt</code></td>
    <td>Amino Acid Sequence String of PDB</td>
  </tr>
  <tr>
    <td><code>{pdb_id}_nowat.pdb</code></td>
    <td>PDB with all water and MG ions removed</td>
  </tr>
  <tr>
    <th colspan="2">🧬 After BIOEMU Inference</th>
  </tr>
  <tr>
    <td><code>{pdb_id}_bioemu.pdb</code></td>
    <td>Final BIOEMU pdb file with multi conformations</td>
  </tr>
  <tr>
    <th colspan="2">🔬 After SAM2 Inference</th>
  </tr>
  <tr>
    <td><code>{pdb_id}_sam2.pdb</code></td>
    <td>Final SAM2 pdb file with multi conformations</td>
  </tr>
  <tr>
    <td><code>sam2.top.pdb</code></td>
    <td>SAM2 Topology PDB File</td>
  </tr>
  <tr>
    <td><code>sam2.traj.dcd</code></td>
    <td>SAM2 Trajectory DCD File</td>
  </tr>
</table>