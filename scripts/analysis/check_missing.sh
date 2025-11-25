# get_missing.sh <dataset_name>
# Returns the Predictors that are missing from the dataset.

#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dataset_name>"
    exit 1
fi

DATASET_NAME=$1
DATASET_TXT="./splits/${DATASET_NAME}.txt"


if [ ! -f "$DATASET_TXT" ]; then
    echo "Dataset file $DATASET_TXT does not exist."
    exit 1
fi

PDB_IDS=$(cat "$DATASET_TXT")

for PDB_ID in $PDB_IDS; do
    oneMissing=0
    if [ ! -f "./PDBs/${PDB_ID}/${PDB_ID}_bioemu.pdb" ]; then
        echo "[$PDB_ID] Missing BioEmu Prediction"
        oneMissing=1
    fi
    if [ ! -f "./PDBs/${PDB_ID}/${PDB_ID}_alphaflow.pdb" ]; then
        echo "[$PDB_ID] Missing AlphaFlow Prediction"
        oneMissing=1
    fi
    if [ ! -f "./PDBs/${PDB_ID}/${PDB_ID}_sam2.pdb" ]; then
        echo "[$PDB_ID] Missing SAM2 Prediction"
        oneMissing=1
    fi
    if [ ! -f "./PDBs/${PDB_ID}/${PDB_ID}_boltz2.pdb" ]; then
        echo "[$PDB_ID] Missing Boltz2 Prediction"
        oneMissing=1
    fi
    if [ ! -f "./PDBs/${PDB_ID}/${PDB_ID}_openfold.pdb" ]; then
        echo "[$PDB_ID] Missing OpenFold Prediction"
        oneMissing=1
    fi

    if [ $oneMissing -eq 0 ]; then

        echo "> [$PDB_ID] All present."
    fi


done