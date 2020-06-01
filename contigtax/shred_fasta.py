#!/usr/bin/env python

import random
from Bio import SeqIO
from argparse import ArgumentParser
import sys


def read_seqs(f):
    return SeqIO.to_dict(SeqIO.parse(f, "fasta"))


def shred(d, prefix=None, existing=False, contigs=10000, minsize=500,
          maxsize=10000):
    """
    Generate random shreds of input fasta file

    :param d: Dictionary of sequences
    :param prefix: Prefix string to append to random contigs
    :param existing: Use existing prefix string ('|' splits prefix)
    :param contigs: Number of contigs to generate
    :param minsize: Minimum size of contigs
    :param maxsize: Maximum size of contigs
    :return: Dictionary of randomly shredded contigs
    """
    random.seed(42)
    shreds = {}
    keys = list(d.keys())

    for i in range(0, contigs):
        # pick a random contig
        key = random.choice(keys)
        if existing:
            prefix = key.split("|")[0]
        if prefix is not None:
            contig_id = ">{}|contig{}".format(prefix, i)
        else:
            contig_id = ">contig{}".format(i)
        keylen = len(d[key])-1
        # pick a random length
        rand_len = random.randrange(minsize, maxsize)
        # if random length is bigger than contig, choose entire contig
        if rand_len >= keylen:
            shreds[contig_id] = str(d[key].seq)
            continue
        # choose whether to start from beginning or end
        if random.choice(["start", "end"]) == "start":
            # if choosing from beginning, pick a random position between
            # the first nucleotide and contig_length - rand_length
            rand_start = random.randrange(0, keylen-rand_len)
            rand_end = rand_start+rand_len
        else:
            rand_end = random.randrange(rand_len, keylen)
            rand_start = rand_end-rand_len
        rand_seq = d[key][rand_start:rand_end]
        shreds[contig_id] = rand_seq.seq
    return shreds


def write_shreds(shreds):
    l = []
    for contig_id in sorted(shreds.keys()):
        seq = shreds[contig_id]
        l.append(len(seq))
        sys.stdout.write("{}\n{}\n".format(contig_id, str(seq)))
    import numpy as np
    sys.stderr.write(
    """
min: {min}
max: {max}
median: {median}
mean: {mean}
    """.format(min=np.min(l), max=np.max(l), median=np.median(l),
               mean=np.mean(l)))


def main(args):
    seqs = read_seqs(args.infile)
    shreds = shred(seqs, args.prefix, args.use_prefix, args.contigs,
                   args.minsize, args.maxsize)
    write_shreds(shreds)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", type=str, help="Input fasta file")
    parser.add_argument("--prefix", type=str, help="prefix to add to ids")
    parser.add_argument("--use-prefix", action="store_true",
                        help="Use already existing prefix for sequences")
    parser.add_argument("--minsize", type=int, default=500,
                        help="Minimum contig size")
    parser.add_argument("--maxsize", type=int, default=10000,
                        help="Maximum contig size")
    parser.add_argument("--contigs", type=int, default=10000,
                        help="Contigs to generate")
    args = parser.parse_args()
    main(args)