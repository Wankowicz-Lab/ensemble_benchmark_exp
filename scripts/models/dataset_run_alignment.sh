#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dataset_name> [max_parallel_jobs]"
    exit 1
fi

DATASET_NAME=$1
MAX_JOBS=${2:-24}

DATASET_TXT="./splits/${DATASET_NAME}.txt"

if [ ! -f "$DATASET_TXT" ]; then
    echo "Dataset $DATASET_TXT does not exist."
    exit 1
fi

wait_for_jobs() {
    while [ $(jobs -p | wc -l) -ge $MAX_JOBS ]; do
        sleep 1
    done
}

mkdir -p "./bin/logs"

for PDBID in $(cat $DATASET_TXT); do
    echo "Processing $PDBID"
    for PREDICTOR in "alphaflow" "sam2" "openfold" "boltz2" "bioemu"; do
        wait_for_jobs 
        
        echo "Starting alignment of $PREDICTOR for $PDBID"
        bash "./scripts/models/complete_phaser_alignment.sh" "$PDBID" "$PREDICTOR" > "./bin/logs/${PDBID}_${PREDICTOR}.log" 2>&1 &
    done
done

echo "Waiting for all alignment jobs to complete..."
wait
echo "All alignment jobs completed!"