from __future__ import print_function
import os,sys
import argparse

import pysam
from collections import defaultdict

import math

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd

import argparse, os


def plot_rand_index_per_error_batch(csv_file, outfile):
    # sns.plt.clf()
    with sns.plotting_context("paper", font_scale=2.4):
        
        rc={'axes.labelsize': 16.0, 'font.size': 16.0, 'legend.fontsize': 10.0, 'axes.titlesize': 14, 'xtick.labelsize': 16.0, 'ytick.labelsize' : 16.0}
        sns.set(rc=rc)
        indata = pd.read_csv(csv_file)
        # ARI,error_rate,nr_samples
        g = sns.pointplot(x="error_rate", y="ARI", hue="dataset", data=indata)

        
        # g = sns.FacetGrid(indata, col="measure_type", size=3, aspect=1.6, col_order=["V-measure", "Completeness", "Homogeneity"], legend_out=True)
        # (g.map(sns.pointplot, "class_size", "measure", "Dataset", hue_order=["ALZ_PB", "BHAM_ONT", "RC0_PB", "HUM_PB", "ZEB_PB", "ENS_PB_100k", "ENS_PB_500k", "ENS_PB_1M"], palette=sns.color_palette("Set1", 8)).despine(left=True).add_legend(title="DATASET", label_order=["ALZ_PB", "BHAM_ONT", "RC0_PB", "HUM_PB", "ZEB_PB", "ENS_PB_100k", "ENS_PB_500k", "ENS_PB_1M"]))
        # g.set(ylim=(0.75, 1))

        # resize figure box to -> put the legend out of the figure
        # box = g.ax.get_position() # get position of figure
        # g.ax.set_position([box.x0, box.y0, box.width * 0.85, box.height]) # resize position

        # # Put a legend to the right side
        # g.ax.legend(loc='center right', bbox_to_anchor=(1.25, 0.5), ncol=1)
        
        # plt.legend(loc='lower right')
        plt.xticks(rotation=90)
        plt.tight_layout()
        # g.set_title('Clustering accuracy per class size')
        # g.set_ylabel("Performance")
        # g.set_xlabel("Class size")
        plt.savefig(args.outfile)
        plt.close()






def main(args):
    plot_rand_index_per_error_batch(args.infile, args.outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--infile', type=str, help='Csv infile with results')
    parser.add_argument('--outfile', type=str, help='Output file with results')
    args = parser.parse_args()

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)

    main(args)




