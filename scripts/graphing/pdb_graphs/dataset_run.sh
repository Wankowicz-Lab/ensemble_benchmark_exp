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


# for prefixing stuff for import order;
get_protein_class() {
    local pdbid=$1
    pdbid=$(echo "$pdbid" | tr '[:upper:]' '[:lower:]')
    
    # CypA proteins
    if [[ "$pdbid" == "4yuo" || "$pdbid" == "3k0m" || "$pdbid" == "5f66" || 
          "$pdbid" == "3k0n" || "$pdbid" == "3k0o" || "$pdbid" == "4yuh" || 
          "$pdbid" == "4yuj" ]]; then
        echo "1-CypA"
    # Lysozyme proteins
    elif [[ "$pdbid" == "2lzt" || "$pdbid" == "6o2h" || "$pdbid" == "4lzt" || 
            "$pdbid" == "2vb1" || "$pdbid" == "1aki" || "$pdbid" == "8dyz" || 
            "$pdbid" == "8dz7" ]]; then
        echo "2-Lysozyme"
    # UBI proteins
    elif [[ "$pdbid" == "5tof" || "$pdbid" == "5tog" ]]; then
        echo "3-UBI"
    # STN proteins
    elif [[ "$pdbid" == "2pyk" || "$pdbid" == "2pzw" ]]; then
        echo "4-STN"
    # All others
    else
        echo "5-other"
    fi
}

for PDBID in $(cat $DATASET_TXT); do
    echo "Running $SCRIPT_NAME for $PDBID..."
    mkdir -p "./bin/graphs/$DATASET_NAME/$PDBID"
    SCRIPT_NAME_WITHOUT_EXT="${SCRIPT_NAME%.*}"
    CLASS_PREFIX=$(get_protein_class "$PDBID")
    
    python "./scripts/graphing/pdb_graphs/$SCRIPT_NAME" "$PDBID" "./bin/graphs/latest/md_bmk/heavy_${CLASS_PREFIX}_$PDBID-$SCRIPT_NAME_WITHOUT_EXT.png"
done

