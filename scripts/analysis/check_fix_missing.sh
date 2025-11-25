# recheck_missing.sh <dataset_name>
# Finds the Predictors that are missing from the dataset and sees if they just didn't get moved to the right place.
# non-moved will likely be failures of phaser alignment.

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
        echo "[$PDB_ID] Missing BioEmu Prediction. Checking if an Ensemble exists..."
        if [ -f "./PDBs/${PDB_ID}/bioemu_bin/${PDB_ID}_ensemble.pdb" ]; then
            echo "[$PDB_ID] BioEmu Ensemble found. Moving to correct location."
            cp "./PDBs/${PDB_ID}/bioemu_bin/${PDB_ID}_ensemble.pdb" "./PDBs/${PDB_ID}/${PDB_ID}_bioemu.pdb"
        else
            echo "[$PDB_ID] No BioEmu Ensemble found."
        fi
        oneMissing=1
    fi
    if [ ! -f "./PDBs/${PDB_ID}/${PDB_ID}_alphaflow.pdb" ]; then
        echo "[$PDB_ID] Missing AlphaFlow Prediction. Checking if an Ensemble exists..."
        if [ -f "./PDBs/${PDB_ID}/alphaflow_bin/${PDB_ID}_ensemble.pdb" ]; then
            echo "[$PDB_ID] AlphaFlow Ensemble found. Moving to correct location."
            cp "./PDBs/${PDB_ID}/alphaflow_bin/${PDB_ID}_ensemble.pdb" "./PDBs/${PDB_ID}/${PDB_ID}_alphaflow.pdb"
        else
            echo "[$PDB_ID] No AlphaFlow Ensemble found."
        fi
        oneMissing=1
    fi
    if [ ! -f "./PDBs/${PDB_ID}/${PDB_ID}_sam2.pdb" ]; then
        echo "[$PDB_ID] Missing SAM2 Prediction. Checking if an Ensemble exists..."
        if [ -f "./PDBs/${PDB_ID}/sam2_bin/${PDB_ID}_ensemble.pdb" ]; then
            echo "[$PDB_ID] SAM2 Ensemble found. Moving to correct location."
            cp "./PDBs/${PDB_ID}/sam2_bin/${PDB_ID}_ensemble.pdb" "./PDBs/${PDB_ID}/${PDB_ID}_sam2.pdb"
        else
            echo "[$PDB_ID] No SAM2 Ensemble found."
        fi
        oneMissing=1
    fi
    if [ ! -f "./PDBs/${PDB_ID}/${PDB_ID}_boltz2.pdb" ]; then
        echo "[$PDB_ID] Missing Boltz2 Prediction. Checking if an Ensemble exists..."
        if [ -f "./PDBs/${PDB_ID}/boltz2_bin/${PDB_ID}_ensemble.pdb" ]; then
            echo "[$PDB_ID] Boltz2 Ensemble found. Moving to correct location."
            cp "./PDBs/${PDB_ID}/boltz2_bin/${PDB_ID}_ensemble.pdb" "./PDBs/${PDB_ID}/${PDB_ID}_boltz2.pdb"
        else
            echo "[$PDB_ID] No Boltz2 Ensemble found."
        fi
        oneMissing=1
    fi
    if [ ! -f "./PDBs/${PDB_ID}/${PDB_ID}_openfold.pdb" ]; then
        echo "[$PDB_ID] Missing OpenFold Prediction. Checking if an Ensemble exists..."
        if [ -f "./PDBs/${PDB_ID}/openfold_bin/${PDB_ID}_ensemble.pdb" ]; then
            echo "[$PDB_ID] OpenFold Ensemble found. Moving to correct location."
            cp "./PDBs/${PDB_ID}/openfold_bin/${PDB_ID}_ensemble.pdb" "./PDBs/${PDB_ID}/${PDB_ID}_openfold.pdb"
        else
            echo "[$PDB_ID] No OpenFold Ensemble found."
        fi 
        oneMissing=1
    fi

    if [ $oneMissing -eq 0 ]; then

        echo "> [$PDB_ID] All present."
    fi


done