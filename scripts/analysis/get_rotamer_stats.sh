#!/bin/bash
# get rotamer stats using phenix.rotalyze

PDB_CODE=$1
PREDICTOR=$2

# if no predictor, we will generate for deposited PDB.
if [ -z "$PDB_CODE" ]; then
    echo "Usage: $0 <pdb_id>"
    exit 1
fi

source /dors/wankowicz_lab/Phenix-2.0/phenix_env.sh
#source /programs/sbgrid.shrc # doesn't work sometimes
set -e

if [ -z "$PREDICTOR" ]; then
    ABS_PDB_FOLDER=$(realpath "./PDBs/${PDB_CODE,,}/")
    ABS_ENSEMBLE_PDB="$ABS_PDB_FOLDER/${PDB_CODE,,}_final.pdb"
    PREDICTOR_BIN="$ABS_PDB_FOLDER"
    ROTALYZE_OUTPUT="$PREDICTOR_BIN/rotalyze"

    mkdir -p $ROTALYZE_OUTPUT
    cd $ROTALYZE_OUTPUT

    phenix.rotalyze model=$ABS_ENSEMBLE_PDB outliers_only=False > out.txt

    cd - 
    echo "Rotamer stats for $PDB_CODE DEPOSITED saved in $ROTALYZE_OUTPUT/out.txt"
else
    ABS_PDB_FOLDER=$(realpath "./PDBs/${PDB_CODE,,}/")
    ABS_ENSEMBLE_PDB="$ABS_PDB_FOLDER/${PDB_CODE,,}_${PREDICTOR}.pdb"
    PREDICTOR_BIN="$ABS_PDB_FOLDER/${PREDICTOR,,}_bin/"
    ROTALYZE_OUTPUT="$PREDICTOR_BIN/rotalyze"

    mkdir -p $ROTALYZE_OUTPUT
    cd $ROTALYZE_OUTPUT

    phenix.rotalyze model=$ABS_ENSEMBLE_PDB outliers_only=False > out.txt

    cd - 
    echo "Rotamer stats for $PDB_CODE with predictor $PREDICTOR saved in $ROTALYZE_OUTPUT/out.txt"
fi
