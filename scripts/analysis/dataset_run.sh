# dataset_run.sh <script_name> <dataset_name> 
# This will loop call the Python script_name for all PDBs in the datase

#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <script_name> <dataset_name>"
    exit 1
fi


SCRIPT_NAME=$1
DATASET_NAME=$2

DATASET_TXT="./splits/${DATASET_NAME}.txt"

if [ ! -f "$DATASET_TXT" ]; then
    echo "Dataset $DATASET_TXT does not exist."
    exit 1
fi

PREDICTOR_TYPES=("bioemu" "sam2" "alphaflow" "boltz2" "openfold")


for PDBID in $(cat $DATASET_TXT); do
    #for PREDICTOR in "${PREDICTOR_TYPES[@]}"; do
    #    echo "Running $SCRIPT_NAME for $PDBID... and $PREDICTOR"
        "./scripts/analysis/$SCRIPT_NAME" "$PDBID" # "$PREDICTOR" #"--mode" "backbone"
    #done
done