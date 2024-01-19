The MPRA IVSM folder contains code for obtaining the IVSM for the FIRE enhancer from our MPRA experiment.

DeepChickenBrain contains a full example of training and using one of our enhancer models.

!!! IMPORTANT FOR CALCULATING CONTRIBUTION SCORES
    Add the following code to $CONDA_PREFIX/lib.python3.8/site-packages/shap/explainers/deep_tf.py at line 280:
        elif output_rank_order.isnumeric():
                    model_output_ranks = np.argsort(-model_output_values)
                    model_output_ranks[0] = int(output_rank_order)

