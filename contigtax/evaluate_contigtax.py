#!/usr/bin/env python

import pandas as pd
from argparse import ArgumentParser
from contigtax.assign import make_lineage_df
import sys


def read_taxfile(f):
    return pd.read_table(f, sep="\t", index_col=0, header=None,
                         names=["seqid", "taxid"])


def evaluate(f, taxmap, ranks):
    """
    Read taxonomy assignment file with taxids and compared to known
    taxonomy

    :param f: contigtax output (--taxidout) with taxids instead of names
    :param taxmap: TSV file mapping query ids to taxids
    :param ranks: Ranks to evaluate
    :return: DataFrame with 0 (incorrect) or 1 (correct) assignment for
             each query and rank
    """
    df = pd.read_csv(f, header=0, sep="\t", index_col=0)
    e = {}
    for rank in ranks:
        # Iterate all assigned ranks
        for q in df.loc[df.loc[:, rank] == taxmap.loc[df.index, rank]].index:
            if q not in e.keys():
                e[q] = {rank: 1}
            else:
                e[q][rank] = 1
        for q in df.loc[df.loc[:, rank] != taxmap.loc[df.index, rank]].index:
            if q not in e.keys():
                e[q] = {rank: 0}
            else:
                e[q][rank] = 0
    return pd.DataFrame(e).T[ranks]


def main():
    parser = ArgumentParser()
    parser.add_argument("infile", type=str,
                        help="contigtax assignment file with taxids (use "
                             "--taxidout with contigtax assign)")
    parser.add_argument("taxfile", type=str,
                        help="File mapping sequence id to true taxonomy id")
    parser.add_argument("-t", "--taxdir", type=str, default="taxonomy",
                        help="Directory with ete3 sqlite database")
    parser.add_argument("--dbname", type=str, default="taxonomy.sqlite",
                        help="Name of sqlite database file")
    parser.add_argument("-r", "--ranks", nargs="+",
                        default=["superkingdom", "phylum", "class", "order",
                                 "family", "genus", "species"])
    args = parser.parse_args()
    taxmap = read_taxfile(args.taxfile)
    taxids = [int(x) for x in taxmap.taxid.unique()]
    lineage_df, name_dict = make_lineage_df(taxids=taxids, taxdir=args.taxdir,
                                            dbname=args.dbname,
                                            ranks=args.ranks)
    taxmap = pd.merge(taxmap, lineage_df, left_on="taxid", right_index=True,
                      how="left")
    eval_df = evaluate(args.infile, taxmap, args.ranks)
    eval_df.to_csv(sys.stdout, sep="\t")
    sys.stderr.write("True positive rate (%):\n")
    (eval_df.sum().div(len(eval_df))*100).to_csv(sys.stderr, sep="\t")


if __name__ == '__main__':
    main()
