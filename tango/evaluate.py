#!/usr/bin/env python
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import pandas as pd
from argparse import ArgumentParser
from tango.assign import make_lineage_df
from multiprocessing import Pool
import tqdm
import sys


def read_taxfile(f):
    return pd.read_table(f, sep="\t", index_col=0, header=None, names=["seqid", "taxid"])


def evaluate(f, taxmap, ranks):
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
                        help="Tango assignment file")
    parser.add_argument("taxfile", type=str,
                        help="File mapping sequence id to true taxonomy id")
    parser.add_argument("-t", "--taxdir", type=str, default="taxonomy",
                        help="Directory with ete3 sqlite database")
    parser.add_argument("--dbname", type=str, default="taxonomy.sqlite",
                        help="Name of sqlite database file")
    parser.add_argument("-r", "--ranks", nargs="+",
                        default=["superkingdom", "phylum", "class", "order", "family", "genus", "species"])
    args = parser.parse_args()
    taxmap = read_taxfile(args.taxfile)
    lineage_df = make_lineage_df(taxmap.taxid.unique(), args.taxdir, args.dbname, args.ranks)
    lineage_df.drop(args.ranks, axis=1, inplace=True)
    lineage_df.rename(columns = lambda x: x.replace(".name",""), inplace=True)
    taxmap = pd.merge(taxmap, lineage_df, left_on="taxid", right_index=True, how="left")
    eval_df = evaluate(args.infile, taxmap, args.ranks)
    eval_df.to_csv(sys.stdout, sep="\t")
    sys.stderr.write("True positive rate (%):\n")
    (eval_df.sum().div(len(eval_df))*100).to_csv(sys.stderr, sep="\t")


if __name__ == '__main__':
    main()
