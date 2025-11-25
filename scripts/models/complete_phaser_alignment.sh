# complete_alignment.sh <pdb_id> <predictor_name>

# Sometimes some issue may cause MDTraj alignment to fail, but to not have to 
# re-run the entire script, we can just run this script to pick up where
# we left off.

# This is for AFTER ./PDBs/{pdb}/[etc]_bin/{pdb}_ensemble.pdb exists,
# but aligned model does not.

set -e 

PDB_ID=$1
PREDICTOR_NAME=$2

if [ -z "$PDB_ID" ] || [ -z "$PREDICTOR_NAME" ]; then
    echo "Usage: $0 <pdb_id> <predictor_name>"
    exit 1
fi

PDB_DIR="./PDBs/${PDB_ID}"
ESB_PATH="${PDB_DIR}/${PREDICTOR_NAME}_bin/${PDB_ID}_ensemble.pdb"
MTZ_PATH="${PDB_DIR}/${PDB_ID}_final.mtz"


if [ ! -f "$ESB_PATH" ]; then
    echo "Error: PDB file does not exist at $ESB_PATH"
    exit 1
fi


TOPOLOGY_PATH="${PDB_DIR}/${PDB_ID}_final.pdb"
python ./scripts/helpers/align_with_mdtraj.py "$ESB_PATH" "$TOPOLOGY_PATH"

bash ./scripts/helpers/align_with_phaser.sh "${PDB_ID,,}" "$PREDICTOR_NAME"

echo "[complete_alignment.sh] Final, aligned PDB file: $PDB_DIR/${PDB_ID,,}_${PREDICTOR_NAME}.pdb"

