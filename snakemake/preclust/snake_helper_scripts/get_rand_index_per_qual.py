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


import networkx as nx
import argparse, os
import pysam

from collections import defaultdict

# from sklearn.metrics.cluster import v_measure_score, completeness_score, homogeneity_score
from sklearn.metrics.cluster import adjusted_rand_score, v_measure_score, completeness_score, homogeneity_score
from sklearn.metrics import fowlkes_mallows_score #, adjusted_mutual_info_score

def parse_inferred_clusters_tsv(tsv_file, args):
    infile = open(tsv_file , "r")
    clusters = {}
    for line in infile:
        cluster_id, read_acc = line.strip().split("\t")
        if args.simulated:
            read_acc = "_".join([item for item in read_acc.split("_")[:-1] ])
        elif args.ont:
            read_acc = read_acc.split("_runid")[0]
        else:
            read_acc = read_acc.split("_strand")[0] 

        clusters[read_acc] = int(cluster_id)
    return clusters


def parse_true_clusters(ref_file):
    classes = defaultdict(dict)
    class_index = 1
    class_ranges = {}
    ref_id_to_chrom = {}
    alignment_counter = defaultdict(int)
    
    prev_chrom = -1
    curr_class_id = -1
    prev_class_start = -1
    prev_class_stop = -1
    prev_read_id = ""
    unique_reads = set()
    unclassified = 0
    for read in ref_file.fetch(until_eof=True):
        unique_reads.add(read.query_name)
        if read.is_unmapped:
            unclassified += 1
            continue
        if read.is_secondary or read.is_supplementary: # deal with supplementary alignments!!
            continue
        # print(read.query_name, read.flag)
        assert prev_read_id != read.query_name

        chrom = read.reference_name
        if chrom != prev_chrom:
            curr_class_id += 1
            classes[read.query_name] = curr_class_id
            prev_chrom = chrom
            prev_class_start = read.reference_start
            prev_class_stop = read.reference_end

        else:
            read_ref_start = read.reference_start
            read_ref_end = read.reference_end
            if read_ref_start > prev_class_stop:
                curr_class_id += 1
                classes[read.query_name] = curr_class_id
                prev_class_start = read.reference_start
                prev_class_stop = read.reference_end
            else:
                classes[read.query_name] = curr_class_id
            

        prev_class_stop = max(read.reference_end, prev_class_stop)
        prev_read_id = read.query_name



        # classes[read.query_name] = int(read.reference_id)  # chrom
        ref_id_to_chrom[int(read.reference_id)] = chrom
        alignment_counter[int(read.reference_id)] += 1
        # if chrom not in class_ranges:
        #     class_ranges[chrom] = {}

        # print(chrom, read_ref_start, read_ref_end)
        # for start, stop in class_ranges[chrom]:
        #     if start <= read_ref_start and  read_ref_end <= stop:
        #         # entirly within
        #     elif

        # classes[read.query_name] =  #read.reference_name.split("|")[0]  #"|".join(read.reference_name.split("|")[:2])
    print(ref_id_to_chrom)
    print()
    print(alignment_counter)
    print()
    return classes, len(unique_reads), unclassified


def parse_true_clusters_simulated(ref_file):
    classes = defaultdict(dict)
    for read in ref_file.fetch(until_eof=True):
        classes[read.query_name] = read.reference_name.split("|")[0] # by gene id 
        # classes[read.query_name] = read.reference_name.split("|")[1] # by transcript id 
    return classes


import numpy as np
from scipy.misc import comb
def rand_index_score(clusters, classes):
    tp_plus_fp = comb(np.bincount(clusters), 2).sum()
    tp_plus_fn = comb(np.bincount(classes), 2).sum()
    A = np.c_[(clusters, classes)]
    tp = sum(comb(np.bincount(A[A[:, 0] == i, 1]), 2).sum()
             for i in set(clusters))
    fp = tp_plus_fp - tp
    fn = tp_plus_fn - tp
    tn = comb(len(A), 2) - tp - fp - fn
    return (tp + tn) / (tp + fp + fn + tn)


def compute_rand_index_per_error_rate(clusters, classes, read_qualities, outfile_path, dataset, tool):
    """
        split reads up into disjoint sets based on error error_rate
        compute the rand indec for each disjoint set of reads. 
        print to tsv file the rand_index, the number of reads in that batch, the error rate of the batch 
        Compute homog, completeness, V measure on each of these subdata sets.
    """

    batch_by_error_rate = { 0.00 : {} ,0.01 : {} ,0.02 : {}, 0.03 : {}, 0.04 : {}, 
                            0.05 : {} ,0.06 : {}, 0.07 : {}, 0.08 : {}, 0.09 : {},
                             0.10 : {} ,0.11 : {}, 0.12 : {}, 0.13 : {},  0.14 : {},
                             0.15 : {} ,0.16 : {}, 0.17 : {}, 0.18 : {}, 0.19 : {},  
                              0.20 : {} ,0.21 : {}, 0.22 : {}, 0.23 : {}, 0.24 : {}, 0.25 : {} }
    print(list(classes.keys())[:10])
    print()
    print(list(clusters.keys())[:10])
    print()
    print(list(read_qualities.keys())[:10])
    cluster_index = max(list(clusters.values())) +1
    # initialize
    for e in batch_by_error_rate:
        batch_by_error_rate[e]["class"] = []
        batch_by_error_rate[e]["cluster"] = []

    for read_acc, error_rate in read_qualities.items():
        if read_acc in classes and read_acc in clusters:
            eps_r = round(error_rate, 2)
            # print(error_rate, eps_r)
            if eps_r > 0.25:
                eps_r = 0.25

            batch_by_error_rate[eps_r]["class"].append(classes[read_acc])
            batch_by_error_rate[eps_r]["cluster"].append(clusters[read_acc])
        
        elif read_acc in classes: # these reads were omitted in clustering output, assign new unique inferred sinleton cluster
            eps_r = round(error_rate, 2)
            # print(error_rate, eps_r)
            if eps_r > 0.25:
                eps_r = 0.25

            batch_by_error_rate[eps_r]["class"].append(classes[read_acc])
            batch_by_error_rate[eps_r]["cluster"].append(cluster_index)
            cluster_index += 1



    outfile = open(outfile_path, "w")
    outfile.write("ARI,RI,FMI,V,C,H,error_rate,nr_samples,dataset,tool\n")
    for e in batch_by_error_rate:
        class_list = batch_by_error_rate[e]["class"]
        cluster_list =  batch_by_error_rate[e]["cluster"]
        if e == 0.0:
            print(class_list)
            print()
            print(cluster_list)
            # take subsample only
            class_list = class_list[:75000]
            cluster_list = cluster_list[:75000]

        ARI = adjusted_rand_score(class_list, cluster_list)
        RI = rand_index_score(class_list, cluster_list)
        FMI = fowlkes_mallows_score(class_list, cluster_list)
        # AMI = adjusted_mutual_info_score(class_list, cluster_list)
        v_score = v_measure_score(class_list, cluster_list)
        compl_score = completeness_score(class_list, cluster_list)
        homog_score = homogeneity_score(class_list, cluster_list)

        nr_samples = len(batch_by_error_rate[e]["class"])
        print(ARI, RI, e, nr_samples, FMI)


        outfile.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}\n".format(ARI, RI, FMI, v_score, compl_score, homog_score, e, nr_samples, dataset, tool))
    outfile.close



    # # bin together into
    # bins= [(0,5), (6,10), (11,20), (21,50), (50, 10000000000)]

    # for  (l, u) in bins:
    #     class_list = []
    #     cluster_list = []
    #     unique_classes = set()
    #     for expr_level in  class_sorted_by_expression:
    #         print(expr_level, l, u)
    #         if l <= expr_level <= u:
    #             for class_id, cluster_id in zip(class_sorted_by_expression[expr_level], cluster_sorted_by_expression[expr_level]):
    #                 class_list.append(class_id)
    #                 cluster_list.append(cluster_id)
    #                 unique_classes.add(class_id)


    #     v_score = v_measure_score(class_list, cluster_list)
    #     compl_score = completeness_score(class_list, cluster_list)
    #     homog_score = homogeneity_score(class_list, cluster_list)
    #     print("Bin size: {0}-{1}:".format(l, u))
    #     print("Nr samples:", len(unique_classes))
    #     print("V:", v_score, "Completeness:", compl_score, "Homogeneity:", homog_score)
    #     if l != 50:
    #         tmp_file_.write("{0},{1},{2},{3}\n".format(str(l)+ "-" + str(u), v_score, len(unique_classes), "V-measure" ))
    #         tmp_file_.write("{0},{1},{2},{3}\n".format(str(l)+ "-" + str(u), compl_score,  len(unique_classes), "Completeness" ))
    #         tmp_file_.write("{0},{1},{2},{3}\n".format(str(l)+ "-" + str(u), homog_score, len(unique_classes), "Homogeneity" ))
    #     else:
    #         tmp_file_.write("{0},{1},{2},{3}\n".format(">" + str(l), v_score, len(unique_classes), "V-measure" ))
    #         tmp_file_.write("{0},{1},{2},{3}\n".format(">" + str(l), compl_score,  len(unique_classes), "Completeness" ))
    #         tmp_file_.write("{0},{1},{2},{3}\n".format(">" + str(l), homog_score, len(unique_classes), "Homogeneity" ))
    # tmp_file_.close

    return outfile.name



'''
    Below code taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].replace(" ", "_"), [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break


def get_read_qualities(read_quals):
    D = {chr(i) : 10**( - (ord(chr(i)) - 33)/10.0 )  for i in range(128)}

    if read_quals[-1] == "q": # is fastq
        read_qualities = {}
        for i, (acc, (seq, qual)) in enumerate(readfq(open(read_quals, 'r'))):
            poisson_mean = sum([ qual.count(char_) * D[char_] for char_ in set(qual)])
            error_rate = poisson_mean/float(len(qual))

            if args.simulated:
                acc = "_".join([item for item in acc.split("_")[:-1] ])
            elif args.ont:
                acc = acc.split("_runid")[0]
            else:
                acc = acc.split("_strand")[0] 

            read_qualities[acc] = error_rate
            # print(acc, error_rate)

    elif  read_quals[-1] == "m": # is bam 
        pass
    return read_qualities


def main(args):

    # read qualities is a dict {read_acc : epsilon_r}
    read_qualities = get_read_qualities(args.read_quals)
    print("done reading read qualities")
    clusters = parse_inferred_clusters_tsv(args.clusters, args)
    if args.simulated:
        ref_file = pysam.AlignmentFile(args.classes, "r", check_sq=False)
        classes  = parse_true_clusters_simulated(ref_file)
        tot_nr_reads = len(classes) # by simulation we know classes of all reads, they are therefore the same number.
    else:
        ref_file = pysam.AlignmentFile(args.classes, "rb", check_sq=False)
        classes, tot_nr_reads, unclassified  = parse_true_clusters(ref_file)

    print("done gettinng classes and clusters")

    compute_rand_index_per_error_rate(clusters, classes, read_qualities, args.outfile, args.dataset, args.tool)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--clusters', type=str, help='Inferred clusters (tsv file)')
    parser.add_argument('--classes', type=str, help='A sorted and indexed bam file.')
    parser.add_argument('--read_quals', type=str, help='A fastq or a bam file.')
    parser.add_argument('--dataset', type=str, help='dataset label name.')
    parser.add_argument('--tool', type=str, help='The tool being evaluated.')
    parser.add_argument('--simulated', action="store_true", help='Simulated data, we can simply read correct classes from the ref field.')
    parser.add_argument('--ont', action="store_true", help='ONT data, parsing accessions differently.')
    parser.add_argument('--outfile', type=str, help='Output file with results')
    args = parser.parse_args()

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)

    main(args)




