#!/bin/bash

# align_with_phaser.sh <pdb_id> <predictor>

#source /dors/wankowicz_lab/Phenix-2.0/phenix_env.sh
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

DEPOSITED_PATH="${PDB_ROOT}/${PDB,,}_final.pdb"
MTZ_PATH="${PDB_ROOT}/${PDB,,}_final.mtz"
ENSEMBLE_PATH="${PR_BIN}/${PDB,,}_ensemble.pdb"
OUTPUT_DIR="${PR_BIN}/phaser_output"

ABSOLUTE_DEPOSITED_PATH=$(realpath "$DEPOSITED_PATH")
ABSOLUTE_ENSEMBLE_PATH=$(realpath "$ENSEMBLE_PATH")
ABSOLUTE_MTZ_PATH=$(realpath "$MTZ_PATH")
ABSOLUTE_OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
ABSOLUTE_BIN_TIMINGS_DIR=$(realpath "./bin/timings")

PHASER_LOG_PATH="${ABSOLUTE_OUTPUT_DIR}/phaser_stdout.log"


#ABSOLUTE_BBOX_PATH=$( realpath "./scripts/helpers/internal/get_bbx.py")

START_MS=$(date +%s%3N)

#DEPOSITED_BOUNDING_BOX=$( python "$ABSOLUTE_BBOX_PATH" "$DEPOSITED_PATH")

#echo "[align_with_phaser.sh] [DEPOSITED]: bounding-box = $DEPOSITED_BOUNDING_BOX"
#X=$(echo $DEPOSITED_BOUNDING_BOX | cut -d',' -f1)
#Y=$(echo $DEPOSITED_BOUNDING_BOX | cut -d',' -f2)
#Z=$(echo $DEPOSITED_BOUNDING_BOX | cut -d',' -f3)
#RANGE=$(echo $DEPOSITED_BOUNDING_BOX | cut -d',' -f4)


mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

phaser << EOF > $PHASER_LOG_PATH 2>&1
TITLe $PDB Ensemble Alignment
MODE MR_AUTO
HKLIn $ABSOLUTE_MTZ_PATH
LABIn F=FP SIGF=SIGFP FREE=FreeR_flag
ENSEmble $PDB PDB $ABSOLUTE_ENSEMBLE_PATH RMS 1.0
SEARch ENSEmble $PDB NUMBER 1 COPIES 1
ROOT pa
EOF


#
#TRANSLATION ORTH
#TRANSLATION VOLUME AROUND
#TRANSLATION POINT $X $Y $Z
#TRANSLATION RANGE 5
NEW_BOUNDING_BOX=$(python "$ABSOLUTE_BBOX_PATH" "pa.1.1.pdb")

NEW_X=$(echo $NEW_BOUNDING_BOX | cut -d',' -f1)
NEW_Y=$(echo $NEW_BOUNDING_BOX | cut -d',' -f2)
NEW_Z=$(echo $NEW_BOUNDING_BOX | cut -d',' -f3)

MAGNITUDE=$(echo "scale=2; sqrt(($X - $NEW_X)^2 + ($Y - $NEW_Y)^2 + ($Z - $NEW_Z)^2)" | bc)
echo "[align_with_phaser.sh] [PHASER]: bounding-box distance = $MAGNITUDE Å"


cd - 

END_MS=$(date +%s%3N)

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



mkdir -p $ABSOLUTE_BIN_TIMINGS_DIR
echo "$PDB,$PREDICTOR,$((END_MS - START_MS))" >> "./bin/timings/align_with_phaser.csv"

echo "[align_with_phaser.sh] Alignment completed successfully at $TARGET_FINAL_PATH"

#ALIGNED_BOUNDING_BOX=$( python "./scripts/helpers/internal/get_bbx.py" "$TARGET_FINAL_PATH")

echo "[align_with_phaser.sh] MDTraj bounding-box: $ALIGNED_BOUNDING_BOX"

#A_X=$(echo $ALIGNED_BOUNDING_BOX | cut -d',' -f1)
#A_Y=$(echo $ALIGNED_BOUNDING_BOX | cut -d',' -f2)
#A_Z=$(echo $ALIGNED_BOUNDING_BOX | cut -d',' -f3)

#MAGNITUDE=$(echo "scale=2; sqrt(($X - $A_X)^2 + ($Y - $A_Y)^2 + ($Z - $A_Z)^2)" | bc)
#echo "[align_with_phaser.sh] Alignment BBX Distance from Deposited magnitude: $MAGNITUDE Å"
