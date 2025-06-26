#!/bin/bash

# align_with_phaser.sh <pdb_id> <predictor>

source /dors/wankowicz_lab/Phenix-2.0/phenix_env.sh
set -e


PDB=$1
PREDICTOR=$2

if [ -z "$PDB" ] || [ -z "$PREDICTOR" ]; then
    echo "[align_with_phaser.sh] Usage: $0 <PDB_ID> <PREDICTOR_NAME>"
    exit 1
fi

PDB_ROOT="./PDBs/${PDB,,}"
PR_BIN="${PDB_ROOT}/${PREDICTOR,,}_bin"

if [ ! -d "$PR_BIN" ]; then
    echo "[align_with_phaser.sh] Error: Predictor directory '$PR_BIN' does not exist."
    exit 1
fi

MTZ_PATH="${PDB_ROOT}/${PDB,,}_final.mtz"
ENSEMBLE_PATH="${PR_BIN}/${PDB,,}_ensemble.pdb"
OUTPUT_DIR="${PR_BIN}/phaser_output"

ABSOLUTE_ENSEMBLE_PATH=$(realpath "$ENSEMBLE_PATH")
ABSOLUTE_MTZ_PATH=$(realpath "$MTZ_PATH")
ABSOLUTE_OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

PHASER_LOG_PATH="${ABSOLUTE_OUTPUT_DIR}/phaser_stdout.log"

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

phaser << EOF > $PHASER_LOG_PATH 2>&1
TITLe $PDB Ensemble Alignment
MODE MR_AUTO
HKLIn $ABSOLUTE_MTZ_PATH
LABIn F=FP SIGF=SIGFP
ENSEmble $PDB PDB $ABSOLUTE_ENSEMBLE_PATH RMS 1.0
SEARch ENSEmble $PDB NUM 1
ROOT pa
EOF

cd - 


ALIGNED_ENSEMBLE_PATH="${ABSOLUTE_OUTPUT_DIR}/pa.1.1.pdb"
ALIGNED_CONFORMATION_PATH="${ABSOLUTE_OUTPUT_DIR}/pa.1.pdb"
TARGET_FINAL_PATH="${PDB_ROOT}/${PDB,,}_${PREDICTOR}.pdb"

if [ ! -f "$ALIGNED_ENSEMBLE_PATH" ] && [ ! -f "$ALIGNED_CONFORMATION_PATH" ]; then
    echo "[align_with_phaser.sh] Error: Aligned ensemble file '$ALIGNED_ENSEMBLE_PATH' does not exist: Phaser was unsuccessful. Phaser Log: $PHASER_LOG_PATH"
    echo "[align_with_phaser.sh] MDTRAJ's alignment will be used in the final result."
    cp "$ENSEMBLE_PATH" "$TARGET_FINAL_PATH"
    exit 0
fi

if [ -f "$ALIGNED_ENSEMBLE_PATH" ]; then
    cp "$ALIGNED_ENSEMBLE_PATH" "$TARGET_FINAL_PATH"
elif [ -f "$ALIGNED_CONFORMATION_PATH" ]; then
    cp "$ALIGNED_CONFORMATION_PATH" "$TARGET_FINAL_PATH"
fi


echo "[align_with_phaser.sh] Alignment completed successfully at $TARGET_FINAL_PATH"