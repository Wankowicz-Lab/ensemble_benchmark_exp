# dataset_run.sh <script_name> <dataset_name> <samples>
# This will loop call the script_name for all PDBs in the dataset, for the scripts that only accept one PDB.

#!/bin/bash

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <script_name> <dataset_name> <samples>"
    exit 1
fi



SCRIPT_NAME=$1
DATASET_NAME=$2
SAMPLES=$3

DATASET_TXT="./splits/${DATASET_NAME}.txt"

if [ ! -f "$DATASET_TXT" ]; then
    echo "Dataset $DATASET_TXT does not exist."
    exit 1
fi


for PDBID in $(cat $DATASET_TXT); do
    echo "Running $SCRIPT_NAME for $PDBID with $SAMPLES samples..."
    bash "./scripts/models/$SCRIPT_NAME" "$PDBID" "$SAMPLES"
done

