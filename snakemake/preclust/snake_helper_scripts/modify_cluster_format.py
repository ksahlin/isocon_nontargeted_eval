import argparse, os

# from collections import defaultdict



def main(args):
    center_to_cluster_mapping = {}
    current_cl_index = 0
    cluster_to_reads = {}
    # reads_to_cluster = {}
    for i, line in enumerate(open(args.input, "r")):
        if i == 0: # header
            continue 
        vals = line.strip().split()
        from_read_acc, to_read_acc = vals[0], vals[1]

        if to_read_acc in center_to_cluster_mapping:
            cluster_id = center_to_cluster_mapping[to_read_acc]
            cluster_to_reads[cluster_id].append(from_read_acc)
        else:
            cluster_id = current_cl_index
            cluster_to_reads[cluster_id] = [from_read_acc]
            center_to_cluster_mapping[to_read_acc] = cluster_id
            current_cl_index += 1

    outfile = open(args.output, "w")
    for cl_id, list_read_acc in cluster_to_reads.items():
        assert len(list_read_acc) >= 1
        for read_acc in list_read_acc:
            outfile.write("{0}\t{1}\n".format(cl_id, read_acc))
    outfile.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('input', type=str, help='Inferred clusters (tsv file)')
    parser.add_argument('output', type=str, help='output cluster format (tsv)')
    args = parser.parse_args()


    main(args)