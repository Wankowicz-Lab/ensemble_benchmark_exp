echo "[download_alphaflow_weights.sh] Downloading AlphaFlow weights..."
mkdir -p ./bin/weights
cd ./bin/weights
wget https://storage.googleapis.com/alphaflow/params/alphaflow_pdb_base_202402.pt
cd -
echo "[download_alphaflow_weights.sh] Download complete."