# DeepBrain

DeepBrain contains an example of using our enhancer models to score and understand any region in any genome for the cell types in our datasets.

If the models or accompanying files are helpful for your research please cite the following publication:

"Deep learning of enhancer codes highlights similarities between mammalian and avian telencephalon cell types"

Nikolai Hecker*, Niklas Kempynck*, David Mauduit, Darina Abaffyová, Ioannis Sarropoulus, Carmen Bravo González-Blas, Sam Dieltiens, Roel Vandepoel, Valerie Christiaens, Elke Leysen, Suresh Poovathingal, Gert Hulselmans, Joris De Wit, Stein Aerts


## Installation

1. Clone the repo:
   ```bash
   git clone https://github.com/aertslab/DeepBrain.git

2. Install libraries
   Create conda environment: 
   ```bash
   conda create -n DeepBrain python=3.8
   conda activate DeepBrain
   ```
   Run the installation script to install the required dependencies and allow for GPU usage + SHAP interpretations.
   ```bash
   chmod +x install.sh
   install.sh
   ```

4. Download the DeepBrain models
   ``` bash
   cd DeepBrain/models
   wget -r -nH -np --cut-dirs=2 https://resources.aertslab.org/papers/DeepBrain/.models/
   ```

