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
        original_indata = pd.read_csv(csv_file)
        indata = original_indata.dropna()
        indata = indata.loc[indata['nr_samples'] > 1000]

        # plot comparing qtclust and CARNAC-LR
        indata = indata.loc[indata['dataset'] == "BHAM_ONT"]
        g = sns.pointplot(x="error_rate", y="FMI", hue="tool", data=indata, palette=sns.color_palette("Set2", 8))
        plt.xticks(rotation=90)
        plt.tight_layout()
        g.set(ylim=(0, 1))
        g.set_title('FMI on BHAM_ONT')
        plt.savefig(args.outfile)
        plt.close()

        # plot only investigating qtclust

        # indata = indata.loc[indata['tool'] == "QUBRIC"]
        # g = sns.pointplot(x="error_rate", y="FMI", hue="dataset", data=indata)
        # plt.xticks(rotation=90)
        # plt.tight_layout()
        # # g.set_title('FMI on BHAM_ONT')
        # plt.savefig(args.outfile)
        # plt.close()


        # indata = indata.loc[indata['dataset'] == "Alzheimer_IsoSeq_2016"]

        # ARI,error_rate,nr_samples
        # g = sns.FacetGrid(indata, col="tool", size=3, aspect=1.6, col_order=["qtclust", "CARNAC-LR"], legend_out=True)
        # (g.map(sns.pointplot, "error_rate", "ARI", "dataset", palette=sns.color_palette("Set1", 8)).despine(left=True).add_legend(title="DATASET"))

        # g = sns.FacetGrid(indata, col="dataset", size=3, aspect=1.6, col_order=["Bham_run1_pass"], legend_out=True)
        # (g.map(sns.pointplot, "error_rate", "FMI", "tool", palette=sns.color_palette("Set1", 8)).despine(left=True).add_legend(title="TOOL"))







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




