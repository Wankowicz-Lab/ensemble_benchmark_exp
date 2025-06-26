#!/bin/bash
# prepare_openfold.sh <split_name>

# note: openfold needs to be installed separately and the path will need to be changed:
set -e
OPENFOLD_PATH="/home/aslamaj/workdir/openfold/bin" # (contains openfold clone)


if [ -z "$1" ]; then
    echo "[prepare_openfold.sh] Usage: $0 <split_name>"
    exit 1
fi


SPLIT_NAME=$1
BIN_ROOT="./bin/openfold_data/$SPLIT_NAME"
ABS_BIN_ROOT=$(realpath $BIN_ROOT)

mkdir -p "./bin/openfold_data"
mkdir -p $BIN_ROOT
rm -rf $BIN_ROOT/*

mkdir -p $BIN_ROOT/inputs
mkdir -p $BIN_ROOT/alignments


# CIF/FASTA Input Preparation
echo "[prepare_openfold.sh] Preparing CIF and FASTA files for $SPLIT_NAME..."

full_split=$(cat "./splits/$SPLIT_NAME.txt" | sort -u)

for PDB_ID in $full_split; do
    echo "[prepare_openfold.sh] Processing PDB ID: $PDB_ID"
    PDB_DIR="./PDBs/$PDB_ID"
 
    # Generate a FASTA   
    SEQ_PATH="$PDB_DIR/${PDB_ID}_seq.txt"
    if [ ! -f "$SEQ_PATH" ]; then
        echo "[prepare_openfold.sh] Error: Sequence file not found for $PDB_ID at $SEQ_PATH"
        continue
    fi
    SEQ=$(cat $SEQ_PATH)
    FASTA_OUT_PATH="$BIN_ROOT/inputs/${PDB_ID}.fasta"
    echo ">${PDB_ID}" > $FASTA_OUT_PATH
    echo "$SEQ" >> $FASTA_OUT_PATH
    echo "[prepare_openfold.sh] FASTA created at: $FASTA_OUT_PATH"

    # Copy CIF
    CIF_PATH="$PDB_DIR/${PDB_ID}_final.cif"
    CIF_OUT_PATH="$BIN_ROOT/inputs/${PDB_ID}.cif"
    if [ -f "$CIF_PATH" ]; then
        cp $CIF_PATH $CIF_OUT_PATH
        echo "[prepare_openfold.sh] CIF copied to: $CIF_OUT_PATH"
    else
        echo "[prepare_openfold.sh] Warning: CIF file not found for $PDB_ID at $CIF_PATH"
    fi


    echo "[prepare_openfold.sh] Completed processing for PDB ID: $PDB_ID"
done



# OpenFold Alignment Generation
echo "[prepare_openfold.sh] Generating alignments with openfold..."

cd $OPENFOLD_PATH || exit 1

RELRESOURCES="./openfold/resources"
export BASE_DATA_DIR=$(realpath $RELRESOURCES)

python ./openfold/scripts/precompute_alignments.py \
    $ABS_BIN_ROOT/inputs \
    $ABS_BIN_ROOT/alignments \
    --uniref90_database_path $BASE_DATA_DIR/uniref90/uniref90.fasta \
    --mgnify_database_path $BASE_DATA_DIR/mgnify/mgy_clusters_2022_05.fa \
    --pdb70_database_path $BASE_DATA_DIR/pdb70/pdb70 \
    --uniclust30_database_path $BASE_DATA_DIR/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    --bfd_database_path $BASE_DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt 

cd - 

echo "[prepare_openfold.sh] Alignment generation completed for $SPLIT_NAME at $BIN_ROOT/alignments"