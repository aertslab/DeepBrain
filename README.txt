## Installation

1. Clone the repo:
   ```bash
   git clone https://github.com/aertslab/DeepBrain/.git


DeepChickenBrain contains an example of using our enhancer models to score and understand any region in any genome for the cell types in our datasets.

!!! IMPORTANT FOR CALCULATING CONTRIBUTION SCORES
    Add the following code to $CONDA_PREFIX/lib.python3.8/site-packages/shap/explainers/deep_tf.py at line 280:
        elif output_rank_order.isnumeric():
                    model_output_ranks = np.argsort(-model_output_values)
                    model_output_ranks[0] = int(output_rank_order)

