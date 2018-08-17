import argparse, os

from collections import defaultdict


def parse_inferred_clusters_tsv(tsv_file, args):
    infile = open(tsv_file , "r")
    clusters = {}
    for line in infile:
        cluster_id, read_acc = line.strip().split("\t")
        if args.simulated:
            read_acc = "_".join([item for item in read_acc.split("_")[:-1] ])
        else:
            read_acc = read_acc.split("_strand")[0] 

        clusters[read_acc] = int(cluster_id)
    return clusters

def main(args):

    clusters = parse_inferred_clusters_tsv(args.clusters, args)
    if args.simulated:
        ref_file = pysam.AlignmentFile(args.classes, "r", check_sq=False)
        classes  = parse_true_clusters_simulated(ref_file)
    else:
        ref_file = pysam.AlignmentFile(args.classes, "rb", check_sq=False)
        classes  = parse_true_clusters(ref_file)

    value = compute_V_measure(clusters, classes)
    get_cluster_information(clusters, classes)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('input', type=str, help='Inferred clusters (tsv file)')
    parser.add_argument('output', type=str, help='output cluster format (tsv)')
    args = parser.parse_args()


    main(args)