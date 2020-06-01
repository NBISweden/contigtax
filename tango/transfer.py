#!/usr/bin/env python

import pandas as pd
from tango.assign import series2df
from multiprocessing import Pool
import tqdm
import sys


def stage_contigs(df):
    """
    Creates a list of results where each item is a DataFrame of ORFs
    belonging to a single contig

    Parameters
    ----------
    df: pandas.DataFrame
        Dataframe with contig ids as index and taxonomy at each rank

    Returns
    -------
    contig_res: list
        List of results to submit for processing
    """

    contigs = sorted(df.index.unique())
    total_contigs = len(contigs)
    contig_res = []
    for contig in tqdm.tqdm(contigs, total=total_contigs, unit=" contigs",
                            ncols=100, desc="Staging contigs"):
        contig_res.append(df.loc[contig])
    return contig_res


def contig_lca(r):
    """

    Parameters
    ----------
    r: pandas.DataFrame
        Results for a single contig

    Returns
    -------
    lca: pandas.DataFrame
        One row of a DataFrame corresponding to a contig
    """

    r = series2df(r)
    contig = list(set(r.index))[0]
    r = r.drop(["id"], axis=1)
    lca = pd.DataFrame(["Unclassified"] * len(r.columns), index=r.columns).T
    lca.index = [contig]
    for rank in [r.columns[x] for x in list(range(len(r.columns) - 1, -1, -1))]:
        rank_taxa = r[rank].unique()
        if len(rank_taxa) == 1:
            return r.loc[r[rank].unique()[0] == rank_taxa]
    return lca


def transfer_taxonomy(df, gff, ignore_unc_rank=None, cpus=1, chunksize=1,
                      orf_df_out=False):
    """
    This function transfers taxonomy from ORFs to contigs by doing an LCA
    on the dataframe.

    Parameters
    ----------
    df: pandas.DataFrame
        Taxonomy dataframe from tango assign
    gff: str
        GFF file or tsv file mapping contigs to ORFs
    ignore_unc_rank: bool
        Should unclassified ORFs be ignored when doing LCA for contigs?

    Returns
    -------
    contig_tax: pandas.DataFrame
        Contig taxonomy assignments from ORFs
    orf_tax: pandas.DataFrame
        ORF taxonomy assignments (transferred back from contigs)
    """

    # Read the gff
    gff_df = pd.read_csv(gff, header=None, sep="\t", comment="#",
                         usecols=[0, 8], names=["contig", "id"])
    # If the last column only contains 'NA' values instead assume that
    # contigs are in 1st column and ORFs in 2nd
    if gff_df.loc[gff_df["id"] == gff_df["id"]].shape[0] > 0:
        ids = ["{}_{}".format(gff_df.loc[i, "contig"],
                              gff_df.loc[i, "id"].split(";")[0].split("_")[-1])
               for i in gff_df.index]
        gff_df.loc[:, "id"] = ids
    else:
        gff_df = pd.read_csv(gff, header=None, sep="\t", usecols=[0, 1],
                             names=["contig", "id"])
    # Merge ORF df with contig map
    merged_df = pd.merge(df, gff_df, left_index=True, right_on="id")
    # Filter out ORFs not classified at minimum rank
    if ignore_unc_rank:
        merged_df = merged_df.loc[merged_df[ignore_unc_rank] != "Unclassified"]
    contigs = sorted(merged_df.contig.unique())
    merged_df = merged_df.set_index("contig")
    total_contigs = len(contigs)
    sys.stderr.write("Transferring "
                     "taxonomy from {} ORFs to "
                     "{} contigs with "
                     "{} cpus\n".format(len(merged_df), total_contigs, cpus))
    if cpus == 1:
        contig_taxa = []
        for contig in tqdm.tqdm(contigs, total=total_contigs, unit=" contigs",
                                ncols=100, desc="Inferring taxonomy"):
            contig_taxa.append(contig_lca(merged_df.loc[contig]))
    else:
        with Pool(processes=cpus) as pool:
            contig_taxa = list(tqdm.tqdm(
                pool.imap(contig_lca, stage_contigs(merged_df),
                          chunksize=chunksize), desc="Inferring taxonomy",
                total=total_contigs, unit=" contigs", ncols=100))
    contig_df = pd.concat(contig_taxa)
    # Transfer taxonomy back to orfs if specified
    if orf_df_out:
        orf_df = pd.merge(contig_df, gff_df, left_index=True, right_on="contig",
                          how="right")
        orf_df = orf_df.set_index("id")
        orf_df = orf_df.drop("contig", axis=1)
        orf_df = orf_df.fillna("Unclassified")
    else:
        orf_df = None
    return contig_df, orf_df
