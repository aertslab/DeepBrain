# DeepBrain

DeepBrain contains an example of using our enhancer models to score and understand any region in any genome for the cell types in our datasets.

## Installation and usage

### 1. Clone the repo
   ```bash
   git clone https://github.com/aertslab/DeepBrain.git
   ```

### 2. Install libraries
   Create and activate conda environment: 
   ```bash
   conda env create -f environment.yml
   conda activate DeepBrain
   ```
   If the installation fails for some reason, another option is to run the steps in the following script for a manual installation of all packages and the environment:
   ```bash
   ./install.sh
   ```
   If you are using a GPU (recommended), and it is not found after installation, a potential fix may be to link an installed libcusolver.so.11 to the correct path:
   ```bash
   #Define CUDA_INSTALL_PATH depending on where it is installed on the local machine
   ln -s $CUDA_INSTALL_PATH/CUDA/11.3.1/lib64/libcusolver.so.11 $(python -c "import tensorflow.python as x; print(x.__path__[0])")/libcusolver.so.10
   ```

### 3. Download the DeepBrain models
   The weights of the models are stored using Git Large File Storage (LFS). To download them, you will need to have installed Git LFS (https://git-lfs.com/). On Linux, you can install Git LFS with the following command if it was not installed yet:
   ```bash
   sudo apt-get install git-lfs 
   ```
   Then the following commands are required after installation to retrieve the model weights:
   ``` bash
   git lfs install
   git lfs pull
   ```
   If Git LFS does not work, you can also download the model weights from Zenodo: https://zenodo.org/records/10868679
### 4. Usage
   Run the notebook DeepBrain_example.ipynb for example usage for predicting on genomic regions, getting contribution scores and calculating correlation between cell types. If you are running JupyterLab, you can make the environment visible by running:
   ```bash
   ipython kernel install --user --name DeepBrain --display-name "DeepBrain"
   ```

## Citation
If the models or accompanying files are helpful for your research please cite the following publication:

[Enhancer-driven cell type comparison reveals similarities between the mammalian and bird pallium](https://www.biorxiv.org/content/10.1101/2024.04.17.589795v1)

Nikolai Hecker*, Niklas Kempynck*, David Mauduit, Darina Abaffyová, Roel Vandepoel, Sam Dieltiens, Ioannis Sarropoulus, Carmen Bravo González-Blas, Elke Leysen, Rani Moors, Gert Hulselmans, Lynette Lim, Joris De Wit, Valerie Christiaens, Suresh Poovathingal, Stein Aerts
