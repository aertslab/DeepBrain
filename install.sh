conda create --name DeepBrain python=3.8
conda activate DeepBrain
pip install tensorflow==2.4.1
conda install -c conda-forge cudatoolkit=11.2 cudnn=8.1.0
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh

conda install -c conda-forge shap==0.39.0
conda install -c bioconda weblogo==3.7.9 pybedtools==0.8.2
conda install matplotlib==3.5.1 seaborn==0.11.2 jupyterlab==4.0.8 numpy==1.19.5 pandas==1.4.3
conda install -c bioconda pybigwig==0.3.18

pip install modisco==0.5.16.0

