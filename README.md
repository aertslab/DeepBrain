# DeepBrain

DeepBrain contains an example of using our enhancer models to score and understand any region in any genome for the cell types in our datasets.

If the models or accompanying files are helpful for your research please cite the following publication:

"TITLE"

Nikolai Hecker, Niklas Kempynck, David Mauduit, Darina Abaffyovà, Ioannis Sarropoulus, Carmen Bravo González-Blas, Sam Dieltiens, Roel Vandepoel, Gert Hulselmans, Suresh Poovathingal, Valerie Christiaens, Stein Aerts

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
   Install the packages listed in install.txt
   
   !!! IMPORTANT FOR CALCULATING CONTRIBUTION SCORES
   Add the following code to
   ```
   $CONDA_PREFIX/lib/python3.8/site-packages/shap/explainers/deep_tf.py
   ```
   at line 280:
   ```
   elif output_rank_order.isnumeric():
      model_output_ranks = np.argsort(-model_output_values)
      model_output_ranks[0] = int(output_rank_order)
   ```

4. Download the DeepBrain models
   ``` bash
   cd DeepBrain/models
   wget -r -nH -np --cut-dirs=2 https://resources.aertslab.org/papers/DeepBrain/.models/
   ```

