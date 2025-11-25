# get_rmsf_baseline.sh <dataset_name>
# Gets RMSF for multiconformer PDB models in the dataset through QFIT's RMSF script.

QFIT_RMSF_PATH="/home/aslamaj/workdir/alphaflow-testing/ensemble-merger-qfit/qfit-3.0/scripts/post/qfit_RMSF.py"
DATASET_NAME=$1

if [ -z "$DATASET_NAME" ]; then
    echo "Usage: $0 <dataset_name>"
    exit 1
fi

DATASET_TXT="./splits/${DATASET_NAME}.txt"

if [ ! -f "$DATASET_TXT" ]; then
    echo "Dataset $DATASET_TXT does not exist."
    exit 1
fi

for PDBID in $(cat $DATASET_TXT); do
    echo "Processing PDBID: $PDBID"
    PDB_DIR="./PDBs/${PDBID,,}"
    TOPOLOGY="$PDB_DIR/${PDBID,,}_final.pdb"
    TOPOLOGY_ABSPATH="$(realpath $TOPOLOGY)"
    cd "$PDB_DIR/analysis" || { echo "Directory $PDB_DIR/analysis does not exist."; continue; }
    python $QFIT_RMSF_PATH "$TOPOLOGY_ABSPATH" "$PDBID"
    cd -
done