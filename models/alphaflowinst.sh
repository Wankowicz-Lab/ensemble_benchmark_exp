conda install nvidia/label/cuda-11.8.0::cuda
conda install nvidia/label/cuda-11.8.0::cuda-cudart-dev
conda install nvidia/label/cuda-11.8.0::libcusparse-dev
conda install nvidia/label/cuda-11.8.0::libcusolver-dev
conda install nvidia/label/cuda-11.8.0::libcublas-dev
ln -s $CONDA_PREFIX/lib/libcudart_static.a $CONDA_PREFIX/lib/libcudart.a