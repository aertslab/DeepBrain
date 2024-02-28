# DeepBrain

DeepBrain contains an example of using our enhancer models to score and understand any region in any genome for the cell types in our datasets.

If the models or accompanying files are helpful for your research please cite the following publication:

"Deep learning of enhancer codes highlights similarities between mammalian and avian telencephalon cell types"

Nikolai Hecker*, Niklas Kempynck*, David Mauduit, Darina Abaffyová, Ioannis Sarropoulus, Carmen Bravo González-Blas, Sam Dieltiens, Roel Vandepoel, Valerie Christiaens, Elke Leysen, Suresh Poovathingal, Gert Hulselmans, Joris De Wit, Stein Aerts


## Installation

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
   If GPUs are not found after installation, a potential fix may be to link an installed libcusolver.so.11 to the correct path:
   ```bash
   #Define CUDA_INSTALL_PATH depending on where it is installed on the local machine
   ln -s $CUDA_INSTALL_PATH/CUDA/11.3.1/lib64/libcusolver.so.11 $(python -c "import tensorflow.python as x; print(x.__path__[0])")/libcusolver.so.10
   ```

### 3. Download the DeepBrain models
   ``` bash
   cd DeepBrain
   wget -r -nH -np --cut-dirs=2 https://resources.aertslab.org/papers/DeepBrain/.models/
   ```
### 4. Usage
   Run the notebook DeepBrain_example.ipynb for example usage for predicting on genomic regions, getting contribution scores and calculating correlation between cell types.

