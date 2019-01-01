
import argparse, os
import pysam
'''
    Below function taken from https://github.com/lh3/readfq/blob/master/readfq.py
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
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break


def parse_inferred_carnac_clusters(carnacfile):
    clusters = {}
    for i, line in enumerate(open(carnacfile, "r")):
        read_indices = line.strip().split()
        clusters[i] = []
        for read_index in read_indices:
            clusters[i].append(int(read_index))

    return clusters

def main(args):
    q_char = chr(33 + args.quality)
    if args.fastq:
        reads = {acc: (seq, q_char*len(seq)) for (acc, seq, qual) in readfq(open(args.fastq, 'r'))}
    elif args.flnc:
        flnc_file = pysam.AlignmentFile(args.flnc, "rb", check_sq=False)
        reads = {}
        for read in flnc_file.fetch(until_eof=True):
            reads[read.qname] = (read.seq, q_char * len(read.seq))


    out_file = open(args.outfile, "w")
    for acc in reads:
        seq, q = reads[acc]
        out_file.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, seq, "+", q))

    out_file.close()

 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align predicted transcripts to transcripts in ensembl reference data base.")
    parser.add_argument('--flnc', type=str, help='input flnc')
    parser.add_argument('--fastq', type=str, help='input fastq')
    parser.add_argument('--quality', type=int, help='A sorted and indexed bam file.')
    parser.add_argument('--outfile', type=str, help='Output file with results')
    args = parser.parse_args()

    # if not os.path.exists(args.outfolder):
    #     os.makedirs(args.outfolder)

    main(args)