#!/usr/bin/env python3
import pandas as pd
import argparse
import os
import matplotlib.pyplot as plt

def make_argument_parser():
    """
    Creates an ArgumentParser to read the options for this script from
    sys.argv
    """
    parser = argparse.ArgumentParser(
        description="Intersect readID-enhancerBC-plasmidBC file with enhancerID-enhancerBC to return enhancerID-plasmidBC assignments.")
    parser.add_argument('--sequencing_sat_dir', '-ssd', type=str, required=True, 
                        help='Path to the folder where sequencing saturation files have been generated')
    parser.add_argument('--regex', '-r', type=str, required=False, default='_sat',
                        help='Expression to identify files containing the saturation values')
    parser.add_argument('--ncol', '-nc', type=int, required=False, default=3,
                        help='Number of columns to use in figure plot')
    parser.add_argument('--nrow', '-nr', type=int, required=False, default=2,
                        help='Number of rows to use in figure plot')
    return parser

def main():
    """
    The main executable function
    """
    parser = make_argument_parser()
    args = parser.parse_args()
    ssd = args.sequencing_sat_dir
    r = args.regex
    nc = args.ncol
    nr = args.nrow
    
    files = [name for name in os.listdir(ssd) if r in name]
    figsize = (6.4*nc, 4.8*nr)
    fig=plt.figure(figsize=figsize)
    i = 1
    
    for file in files:
        plt.subplot(nr, nc, i)
        toplot=pd.read_csv(os.path.join(ssd,file), sep=' ', engine='python')
        plt.plot(toplot.iloc[:,0], toplot.iloc[:,1])
        plt.xlabel("Number of fragments")
        plt.ylabel("Unique fragments")
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0), useOffset=True)
        plt.title(file, fontsize=10,  y=1.08)
        i += 1

    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    fig.savefig(os.path.join(ssd,'saturation_plot.pdf'))
    plt.show()
if __name__ == "__main__":
    main()
