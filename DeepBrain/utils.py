import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import tensorflow as tf
import shap

def load_model(output_dir,model_name):
    """
    Loads a trained model from the specified output directory.

    Args:
        output_dir (str): The directory where the model files are stored.
        name_hdf5 (str): The name of the HDF5 file containing the model weights.

    Returns:
        model (tf.keras.Model): The loaded model.
    """
    model_json_file = open(os.path.join(output_dir, model_name+'.json'))
    model_json = model_json_file.read()
    model = tf.keras.models.model_from_json(model_json)
    model.load_weights(os.path.join(output_dir, model_name+'.hdf5'))
    return model

def one_hot_encode_along_row_axis(sequence):
    """
    Converts a DNA sequence to a one-hot encoded array along the row axis.

    Args:
        sequence (str): The DNA sequence to be converted.

    Returns:
        to_return (numpy.ndarray): The one-hot encoded array.
    """
    to_return = np.zeros((1, len(sequence), 4), dtype=np.int8)
    seq_to_one_hot_fill_in_array(zeros_array=to_return[0], sequence=sequence, one_hot_axis=1)
    return to_return


def seq_to_one_hot_fill_in_array(zeros_array, sequence, one_hot_axis):
    """
    Converts a DNA sequence to a one-hot encoded array and fills in the provided zeros_array.

    Args:
        zeros_array (numpy.ndarray): The array to be filled with the one-hot encoded sequence.
        sequence (str): The DNA sequence to be converted.
        one_hot_axis (int): The axis along which the one-hot encoding should be applied. 
                            Must be either 0 or 1.

    Raises:
        AssertionError: If the one_hot_axis is not 0 or 1, or if the shape of the zeros_array 
                        does not match the length of the sequence along the specified axis.

    Returns:
        None
    """
    assert one_hot_axis == 0 or one_hot_axis == 1

    if one_hot_axis == 0:
        assert zeros_array.shape[1] == len(sequence)
    elif one_hot_axis == 1:
        assert zeros_array.shape[0] == len(sequence)

    char_mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': None, 'a': 0, 'c': 1, 'g': 2, 't': 3, 'n': None}

    for i, char in enumerate(sequence):
        char_idx = char_mapping.get(char)
        if char_idx is None:
            continue
        if one_hot_axis == 0:
            zeros_array[char_idx, i] = 1
        elif one_hot_axis == 1:
            zeros_array[i, char_idx] = 1

def plot_model_predictions(model, location, seq_onehot, model_classes):
    """
    Plots the model predictions for a given location.

    Args:
        model (keras.Model): The trained model.
        location (str): The location for which predictions are made.
        seq_onehot (numpy.ndarray): The one-hot encoded sequence data.
        model_classes (list): The list of model classes.

    Returns:
        numpy.ndarray: The model predictions.
    """

    location = location[0]+':'+str(location[1])+'-'+str(location[2])
    prediction = model.predict(seq_onehot)[0]
    plt.figure(figsize=(len(model_classes)*1.3,2))
    plt.ylim([0,1])
    plt.ylabel("model prediction")
    plt.xlabel("cell type")
    plt.title("Model predictions for " + str(location))
    plt.bar(model_classes, prediction)
    plt.show()
    return prediction

def reverse_complement_one_hot_sequence(sequence):
    """
    Reverse complements a one-hot encoded sequence.

    Args:
        sequence (numpy.ndarray): The one-hot encoded sequence.

    Returns:
        reverse_complement (numpy.ndarray): The reverse complement of the input sequence.
    """
    reverse_complement = sequence[::-1, ::-1]
    return reverse_complement

def plot_deepexplainer_givenax(explainer, fig, ntrack, track_no, seq_onehot, cell_type, class_names, plot=True):
    """
    Plots the DeepExplainer weights for a given cell type on a given axis.

    Parameters:
    explainer (DeepExplainer): The DeepExplainer object.
    fig (matplotlib.figure.Figure): The figure object to plot on.
    ntrack (int): The number of tracks in the figure.
    track_no (int): The track number to plot on.
    seq_onehot (numpy.ndarray): The one-hot encoded sequence.
    cell_type (str): The cell type for which the weights are plotted.
    class_names (list): The list of class names.

    Returns:
    matplotlib.axes.Axes: The axes object containing the plot.
    numpy.ndarray: The array of SHAP values.
    """
    target_class = class_names.index(cell_type)
    shap_values_, indexes_ = shap_values(explainer.explainer, seq_onehot,
                                                   output_rank_order=str(target_class),
                                                   ranked_outputs=1,
                                                   check_additivity=False)
    if plot:
        _, ax1 = plot_weights(shap_values_[0][0]*seq_onehot,
                            fig, ntrack, 1, track_no,
                            title='Region explanation for ' + cell_type, subticks_frequency=20, ylab="DeepExplainer")
    else:
        ax1 = None
    shaps = shap_values_[0]*seq_onehot
    shaps = shaps[np.where(shaps!=0)]
    return ax1, shaps


# Adapted from https://github.com/shap/shap/blob/master/shap/explainers/_deep/deep_tf.py
def shap_values(explainer, X, ranked_outputs=None, output_rank_order="max", check_additivity=True):
        # check if we have multiple inputs
        if not explainer.multi_input:
            if type(X) == list and len(X) != 1:
                assert False, "Expected a single tensor as model input!"
            elif type(X) != list:
                X = [X]
        else:
            assert type(X) == list, "Expected a list of model inputs!"
        assert len(explainer.model_inputs) == len(X), "Number of model inputs (%d) does not match the number given (%d)!" % (len(explainer.model_inputs), len(X))

        # rank and determine the model outputs that we will explain
        if ranked_outputs is not None and explainer.multi_output:
            if not tf.executing_eagerly():
                model_output_values = explainer.run(explainer.model_output, explainer.model_inputs, X)
            else:
                model_output_values = explainer.model(X)

            if output_rank_order == "max":
                model_output_ranks = np.argsort(-model_output_values)
            elif output_rank_order == "min":
                model_output_ranks = np.argsort(model_output_values)
            elif output_rank_order == "max_abs":
                model_output_ranks = np.argsort(np.abs(model_output_values))
            elif output_rank_order.isnumeric(): ### Added from the original version
                model_output_ranks = np.argsort(-model_output_values)
                model_output_ranks[0] = int(output_rank_order)
            else:
                assert False, "output_rank_order must be max, min, or max_abs!"
            model_output_ranks = model_output_ranks[:,:ranked_outputs]
        else:
            model_output_ranks = np.tile(np.arange(len(explainer.phi_symbolics)), (X[0].shape[0], 1))

        # compute the attributions
        output_phis = []
        for i in range(model_output_ranks.shape[1]):
            phis = []
            for k in range(len(X)):
                phis.append(np.zeros(X[k].shape))
            for j in range(X[0].shape[0]):
                if (hasattr(explainer.data, '__call__')):
                    bg_data = explainer.data([X[l][j] for l in range(len(X))])
                    if type(bg_data) != list:
                        bg_data = [bg_data]
                else:
                    bg_data = explainer.data

                # tile the inputs to line up with the background data samples
                tiled_X = [np.tile(X[l][j:j+1], (bg_data[l].shape[0],) + tuple([1 for k in range(len(X[l].shape)-1)])) for l in range(len(X))]

                # we use the first sample for the current sample and the rest for the references
                joint_input = [np.concatenate([tiled_X[l], bg_data[l]], 0) for l in range(len(X))]

                # run attribution computation graph
                feature_ind = model_output_ranks[j,i]
                sample_phis = explainer.run(explainer.phi_symbolic(feature_ind), explainer.model_inputs, joint_input)

                # assign the attributions to the right part of the output arrays
                for l in range(len(X)):
                    phis[l][j] = (sample_phis[l][bg_data[l].shape[0]:] * (X[l][j] - bg_data[l])).mean(0)

            output_phis.append(phis[0] if not explainer.multi_input else phis)

        # check that the SHAP values sum up to the model output
        if check_additivity:
            if not tf.executing_eagerly():
                model_output = explainer.run(explainer.model_output, explainer.model_inputs, X)
            else:
                model_output = explainer.model(X)
            for l in range(len(explainer.expected_value)):
                if not explainer.multi_input:
                    diffs = model_output[:, l] - explainer.expected_value[l] - output_phis[l].sum(axis=tuple(range(1, output_phis[l].ndim)))
                else:
                    diffs = model_output[:, l] - explainer.expected_value[l]
                    for i in range(len(output_phis[l])):
                        diffs -= output_phis[l][i].sum(axis=tuple(range(1, output_phis[l][i].ndim)))
                assert np.abs(diffs).max() < 1e-2, "The SHAP explanations do not sum up to the model's output! This is either because of a " \
                                                   "rounding error or because an operator in your computation graph was not fully supported. If " \
                                                   "the sum difference of %f is significant compared the scale of your model outputs please post " \
                                                   "as a github issue, with a reproducable example if possible so we can debug it." % np.abs(diffs).max()

        if not explainer.multi_output:
            return output_phis[0]
        elif ranked_outputs is not None:
            return output_phis, model_output_ranks
        else:
            return output_phis


## Plotting functions from DeepLIFT: https://github.com/kundajelab/deeplift/

def plot_a(ax, base, left_edge, height, color):
    a_polygon_coords = [
        np.array([
            [0.0, 0.0],
            [0.5, 1.0],
            [0.5, 0.8],
            [0.2, 0.0],
        ]),
        np.array([
            [1.0, 0.0],
            [0.5, 1.0],
            [0.5, 0.8],
            [0.8, 0.0],
        ]),
        np.array([
            [0.225, 0.45],
            [0.775, 0.45],
            [0.85, 0.3],
            [0.15, 0.3],
        ])
    ]
    for polygon_coords in a_polygon_coords:
        ax.add_patch(matplotlib.patches.Polygon((np.array([1, height])[None, :] * polygon_coords
                                                 + np.array([left_edge, base])[None, :]),
                                                facecolor=color, edgecolor=color))


def plot_c(ax, base, left_edge, height, color):
    ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=1.3, height=height,
                                            facecolor=color, edgecolor=color))
    ax.add_patch(
        matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=0.7 * 1.3, height=0.7 * height,
                                   facecolor='white', edgecolor='white'))
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 1, base], width=1.0, height=height,
                                              facecolor='white', edgecolor='white', fill=True))


def plot_g(ax, base, left_edge, height, color):
    ax.add_patch(matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=1.3, height=height,
                                            facecolor=color, edgecolor=color))
    ax.add_patch(
        matplotlib.patches.Ellipse(xy=[left_edge + 0.65, base + 0.5 * height], width=0.7 * 1.3, height=0.7 * height,
                                   facecolor='white', edgecolor='white'))
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 1, base], width=1.0, height=height,
                                              facecolor='white', edgecolor='white', fill=True))
    ax.add_patch(
        matplotlib.patches.Rectangle(xy=[left_edge + 0.825, base + 0.085 * height], width=0.174, height=0.415 * height,
                                     facecolor=color, edgecolor=color, fill=True))
    ax.add_patch(
        matplotlib.patches.Rectangle(xy=[left_edge + 0.625, base + 0.35 * height], width=0.374, height=0.15 * height,
                                     facecolor=color, edgecolor=color, fill=True))


def plot_t(ax, base, left_edge, height, color):
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge + 0.4, base],
                                              width=0.2, height=height, facecolor=color, edgecolor=color, fill=True))
    ax.add_patch(matplotlib.patches.Rectangle(xy=[left_edge, base + 0.8 * height],
                                              width=1.0, height=0.2 * height, facecolor=color, edgecolor=color,
                                              fill=True))

default_colors = {0: 'green', 1: 'blue', 2: 'orange', 3: 'red'}
default_plot_funcs = {0: plot_a, 1: plot_c, 2: plot_g, 3: plot_t}

def plot_weights_given_ax(ax, array, height_padding_factor, length_padding, subticks_frequency, highlight, colors=default_colors, plot_funcs=default_plot_funcs):
    if len(array.shape) == 3:
        array = np.squeeze(array)
    assert len(array.shape) == 2, array.shape
    if array.shape[0] == 4 and array.shape[1] != 4:
        array = array.transpose(1, 0)
    assert array.shape[1] == 4
    
    max_pos_height = 0.0
    min_neg_height = 0.0
    heights_at_positions = []
    depths_at_positions = []
    
    for i in range(array.shape[0]):
        acgt_vals = sorted(enumerate(array[i, :]), key=lambda x: abs(x[1]))
        positive_height_so_far = 0.0
        negative_height_so_far = 0.0
        
        for letter in acgt_vals:
            plot_func = plot_funcs[letter[0]]
            color = colors[letter[0]]
            
            if letter[1] > 0:
                height_so_far = positive_height_so_far
                positive_height_so_far += letter[1]
            else:
                height_so_far = negative_height_so_far
                negative_height_so_far += letter[1]
            
            plot_func(ax=ax, base=height_so_far, left_edge=i, height=letter[1], color=color)
        
        max_pos_height = max(max_pos_height, positive_height_so_far)
        min_neg_height = min(min_neg_height, negative_height_so_far)
        heights_at_positions.append(positive_height_so_far)
        depths_at_positions.append(negative_height_so_far)
    
    for color in highlight:
        for start_pos, end_pos in highlight[color]:
            assert start_pos >= 0.0 and end_pos <= array.shape[0]
            min_depth = np.min(depths_at_positions[start_pos:end_pos])
            max_height = np.max(heights_at_positions[start_pos:end_pos])
            ax.add_patch(
                matplotlib.patches.Rectangle(xy=[start_pos, min_depth],
                                             width=end_pos - start_pos,
                                             height=max_height - min_depth,
                                             edgecolor=color, fill=False))
    
    ax.set_xlim(-length_padding, array.shape[0] + length_padding)
    ax.xaxis.set_ticks(np.arange(0.0, array.shape[0] + 1, subticks_frequency))
    height_padding = max(abs(min_neg_height) * (height_padding_factor), abs(max_pos_height) * (height_padding_factor))
    ax.set_ylim(min_neg_height - height_padding, max_pos_height + height_padding)
    return ax


def plot_weights(array, fig, n, n1, n2, title='', ylab='', height_padding_factor=0.2,length_padding=1.0, subticks_frequency=20, colors=default_colors, plot_funcs=default_plot_funcs, highlight={}):
    ax = fig.add_subplot(n, n1, n2)
    ax.set_title(title)
    ax.set_ylabel(ylab)
    y = plot_weights_given_ax(ax=ax, array=array,
                              height_padding_factor=height_padding_factor,
                              length_padding=length_padding,
                              subticks_frequency=subticks_frequency,
                              colors=colors,
                              plot_funcs=plot_funcs,
                              highlight=highlight)
    return fig, ax

## end plotting functions from DeepLIFT: https://github.com/kundajelab/deeplift/


