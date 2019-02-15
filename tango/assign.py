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
import numpy as np
from multiprocessing import Pool
import tqdm
import sys
import gzip as gz
from tango.prepare import init_sqlite_taxdb


def calculate_norm_id(df, nolen=False):
    """
    Calculates normalized identity for diamond results
    :param df: pandas DataFrame of diamond results
    :param nolen: boolean specifying if percent id of a hit is to be normalized by length
    :return: pandas DataFrame with an extra column 'normid'

    If there is a column specifying length of either the query or the subject and if the user has not
    explicitly set 'nolen' to True then the percent identity is normalized by the fraction of aligned positions:

                                     alignment length
            normid = percent id x -----------------------
                                  subject or query length
    If nolen = True then normid = percent id / 100
    """
    # Normalize identity to aligned fraction if subject length is available
    if "slen" in df.columns and nolen == False:
        # Create a new column that is: <alignment length>/<subject length> * <alignment identity>/100
        df = df.assign(lfrac=df.length.div(df.slen))
        # For hits where alignment length is greater than subject length, take the inverse
        df.loc[df.lfrac > 1, "lfrac"] = 1 / df.loc[df.lfrac > 1, "lfrac"]
        df = df.assign(normid=df.lfrac.multiply(df.pident.div(100)))
        df.drop("lfrac", axis=1, inplace=True)
    # If not, just return identity
    else:
        df = df.assign(normid=df.pident.div(100))
    return df


def get_thresholds(df, top=10):
    """
    Here bit-score thresholds are calculated per query an returned in a dictionary
    :param df: pandas DataFrame of diamond results
    :param top: Percentage range of top bitscore
    :return: dictionary with queries as keys and bitscore thresholds as values

    The pandas DataFrame is first sorted by bitscore (high to low), then grouped by query, then for the first entry
    per query the top% of the best hit is calculated and converted to dictionary.
    """
    thresholds = (df.sort_values("bitscore", ascending=False).groupby(level=0).first().bitscore * (
        (100 - top)) / 100).to_dict()
    return thresholds


def get_rank_thresholds(ranks, thresholds):
    """
    Constructs dictionary of rank-specific thresholds

    :param ranks: Taxonomic ranks to assign
    :param thresholds: Thresholds for taxonomic ranks
    :return: Dictionary of thresholds
    """
    rank_thresholds = {}
    assert (len(thresholds) == len(ranks))
    for i, rank in enumerate(ranks):
        rank_thresholds[rank] = float(thresholds[i])
    return rank_thresholds


def add_names(x, taxid, ncbi_taxa):
    """
    This function translates taxonomy ids to names. It operates per-row in the lineage dataframe.

    :param x: DataFrame of one taxid and its taxonomic ranks
    :param taxid: Taxid being evaluated
    :param ncbi_taxa: The ete3 sqlite database connection
    :return: The original DataFrame merged with the taxa names
    """
    # Get a names dictionary for all taxids in the row
    names = ncbi_taxa.get_taxid_translator(list(x.loc[taxid].values) + [taxid])
    n = {}
    # Iterate ranks
    for rank in list(x.columns):
        # Get taxid for the current rank
        t = x.loc[taxid, rank]
        # If taxid is negative it means that there is now classified taxonomy at this rank
        # Instead we get the last known name in the hierarchy and use that name with the 'Unclassified' prefix
        # unless the name is 'root' in which case we just use 'Unclassified'
        if t < 0:
            known_name = names[-t]
            if known_name == "root":
                name = "Unclassified"
            else:
                name = "{}.{}".format("Unclassified", known_name)
        # If taxid is positive we just use the name from the dictionary
        else:
            name = names[t]
        # Add name to a dictionary with keys in the form of {rank}.name
        n["{}.name".format(rank)] = name
    name_df = pd.DataFrame(n, index=[taxid])
    return pd.merge(x, name_df, left_index=True, right_index=True)


def propagate_lower(x, taxid, ranks):
    """
    Shift known ranks down through the taxonomic hierarchy.

    :param x: DataFrame of one taxid and its taxonomic ranks
    :param taxid: Taxid being evaluated
    :param ranks: Ranks used for assigning
    :return: pandas DataFrame updated with missing ranks

    Some proteins in the database may map to a taxonomic rank above the lowest taxonomic rank that we are trying to
    assign. For instance, if we use the ranks 'superkingdom phylum genus species' and a protein maps to a taxid at
    rank phylum then we want to add the taxonomic information at the genus and species levels. This is done here by
    adding the negative taxid of the lowest known rank to the lower ranks.

    Example:
    In the Uniref90 database the entry 'E1GVX1' maps to taxonomy id 838 (rank: genus, name: Prevotella).
    When creating the lineage for taxid 838 we add '-838' to rank species.
    """
    rev_ranks = [ranks[x] for x in list(range(len(ranks) - 1, -1, -1))]
    missing = {}
    known = taxid
    for rank in rev_ranks[0:]:
        if rank not in x.columns:
            missing[rank] = -known
        else:
            known = x.loc[taxid, rank]
    return pd.merge(x, pd.DataFrame(missing, index=[taxid]), left_index=True, right_index=True)


def get_lca(r, ranks):
    """
    Assign lowest common ancestor from a set of taxids.

    :param r: Results for a single query, extracted from the main diamond results file
    :param ranks: Taxonomic ranks to assign taxonomy for
    :return: a tuple of dictionaries with ranks as keys and taxa names/ids as values

    This function takes a query-slice of the diamond results after filtering by score (and rank-threshold if tango mode
    is 'rank_lca' or 'rank_vote'). It then iterates through each rank in reverse order checks how many unique taxids are
    found at that rank. If there's only one taxid
    """
    query = r.index.unique()[0]
    rev_ranks = [ranks[x] for x in list(range(len(ranks) - 1, -1, -1))]
    for i, rank in enumerate(rev_ranks):
        higher_ranks = rev_ranks[i:]
        higher_rank_names = ["{}.name".format(x) for x in higher_ranks]
        c = r.groupby(rank).count()
        if len(c) == 1:
            if len(r) == 1:
                lca = r.loc[query, higher_rank_names].values
                lca_taxids = r.loc[query, higher_ranks].values
            else:
                lca = r.loc[query, higher_rank_names].values[0]
                lca_taxids = r.loc[query, higher_ranks].values[0]
            return dict(zip(higher_ranks, lca)), dict(zip(higher_ranks, lca_taxids))
    return {}, {}


def parse_with_rank_thresholds(r, ranks, rank_thresholds, mode, vote_threshold):
    """Performs parsing by iterating through the ranks in reverse, attempting to assign a taxonomy to lowest rank"""
    # Start from lowest rank
    rev_ranks = [ranks[x] for x in list(range(len(ranks) - 1, -1, -1))]
    for i, rank in enumerate(rev_ranks, start=0):
        # Make sure that LCA is not set below current rank
        allowed_ranks = rev_ranks[i:]
        # Get rank threshold
        threshold = rank_thresholds[rank]
        # Filter results by rank threshold
        try:
            _r = r.loc[r.normid >= threshold]
        except KeyError:
            continue
        if len(_r) == 0:
            continue
        lca = {}
        lca_taxids = {}
        # After filtering, either calculate lca from all filtered taxids

        if mode == "rank_lca":
            lca, lca_taxids = get_lca(_r, allowed_ranks)
        # Or at each rank, get most common taxid
        elif mode == "rank_vote":
            vote = get_rank_vote(_r, rank, vote_threshold)
            if len(vote) > 0:
                lca, lca_taxids = get_lca(vote, allowed_ranks)
        if len(lca.keys()) > 0:
            return lca, lca_taxids
    return {}, {}


def get_rank_vote(r, rank, vote_threshold=0.5):
    """Counts unique taxid from filtered dataframe and sums to current rank"""
    # Create dataframe for unique taxids filtered at this rank threshold
    taxid_counts = pd.DataFrame(dict.fromkeys(r.staxids.unique(), 1), index=["count"]).T
    # Add taxid for rank being investigated
    rank_df = r.groupby("staxids").first().reset_index()[[rank, "staxids"]].set_index("staxids")
    rank_df = pd.merge(taxid_counts, rank_df, left_index=True, right_index=True)
    # Sum counts for current rank
    rank_sum = rank_df.groupby(rank).sum()
    rank_norm = rank_sum.div(rank_sum.sum())
    rank_norm = rank_norm.sort_values("count", ascending=False)
    votes = rank_norm.loc[rank_norm["count"] > vote_threshold]
    if len(votes) > 0:
        return r.loc[r[rank].isin(votes.index)]
    return []


def propagate_taxa(res, ranks):
    """Sets unclassified lower ranks from best assignment at higher ranks"""
    known = ""
    for rank in ranks:
        if res[rank] == "Unclassified":
            if known != "" and "Unclassified" not in known:
                res[rank] = "{}.{}".format("Unclassified", known)
            elif known != "" and "Unclassified" in known:
                res[rank] = known
            else:
                continue
        else:
            known = res[rank]
    return res


def series2df(df):
    """Converts pandas series to pandas dataframe"""
    if str(type(df)) == "<class 'pandas.core.series.Series'>":
        df = pd.DataFrame(df).T
    return df


def read_lengthmap(args, res):
    """
    Reads protein id to protein length map file

    :param args: Input arguments
    :param res: Diamond results dictionary
    :return: Dictionary of protein id to protein length
    """
    open_function = open
    # If lengths are for queries we extract all query ids from the results
    if args.querylenmap:
        mapfile = args.querylenmap
        ids = res.keys()
    # If lengths are for subjects we extract all subject ids from the results
    elif args.subjectlenmap:
        mapfile = args.subjectlenmap
        ids = list(set(res[q][i][0] for q in res.keys() for i in range(0, len(res[q]))))
    if ".gz" in mapfile:
        open_function = gz.open
    # Set lengthmap dictionary here
    lengthmap = dict.fromkeys(ids, 0)
    with open_function(mapfile, 'rt') as fhin:
        for line in tqdm.tqdm(fhin, desc="Reading lengths from {}".format(mapfile), ncols=100, unit=" lines"):
            protid, length = (line.rstrip()).rsplit()
            # Attempt to add length to dictionary for protein and subtract the existing length in case the same
            # id is listed more than once in the file
            # If we get a KeyError this is not the id we are looking for
            # If we get a ValueError then it's not a valid length
            try:
                lengthmap[protid] += int(length) - lengthmap[protid]
            except KeyError:
                continue
            except ValueError:
                continue
    return pd.DataFrame(lengthmap, index=["slen"]).T


def read_taxidmap(f, ids):
    """
    Reads the protein to taxid map file and stores mappings

    :param f: Input file
    :param ids: Unique protein ids to store taxids for
    :return: Dictionary of protein ids to taxid and all unique taxids
    """
    taxidmap = dict.fromkeys(ids, -1)
    open_function = open
    if ".gz" in f:
        open_function = gz.open
    with open_function(f, 'rt') as fhin:
        for line in tqdm.tqdm(fhin, desc="Reading idmap {}".format(f), ncols=100, unit=" lines"):
            items = (line.rstrip()).rsplit()
            # If file has only two columns, assume taxid in second
            if len(items) == 2:
                protid, taxid = items
            # Otherwise, assume format is same as NCBI protein mapping
            else:
                protid, taxid = items[0], items[2]
            # Add map to dictionary
            # We initialize the dictionary with -1 so we make an attempt to add the taxid + 1
            # If the protid is not in the dictionary we skip it
            try:
                taxidmap[protid] += int(taxid) + 1
            except KeyError:
                continue
            except ValueError:
                continue
    return pd.DataFrame(taxidmap, index=["staxids"]).T, list(set(taxidmap.values()))


def read_df(file, taxidmap):
    """
    Reads the blast output into a format that can be processed with the multiprocessing module

    If the diamond search has handled internally then the output looks like this
    query1     subject1   93.6    47      3       0       146     6       79      125     8.5e-16 91.3    314295  128
    query1     subject2  100.0   44      0       0       137     6       484     527     2.5e-15 89.7    9347    530
    query2     subject3      53.5    241     84      2       645     7       15      255     1.3e-53 216.9   864142  279

    where the last two columns are taxid and subject length, respectively.

    Otherwise the output may have the typical blast format 6 output.
    """
    open_function = open
    if ".gz" in file:
        open_function = gz.open
    r = {}
    # f = 0 (tango) or 1 (blast) input
    f = 0
    taxids = []
    with open_function(file, 'rt') as fhin:
        for line in tqdm.tqdm(fhin, desc="Reading {}".format(f), ncols=100, unit=" lines"):
            line = line.rstrip()
            items = line.rsplit()
            # If file has 14 columns we assume it's internal tango standard format
            # with taxid in second to last column and subject length in last column
            if len(items) == 14:
                taxids.append(items[-2])
                f = 0
            # If file has 12 columns we assume it's typical blast 6 format
            # In that case a protein -> taxonomy id map has to be present
            elif len(items) == 12:
                if not taxidmap:
                    sys.exit(
                        "ERROR: Standard blast input detected with no protein -> taxid file specified (--taxidmap).")
                f = 1
            else:
                sys.exit("ERROR: Unrecognized input format")
            try:
                r[items[0]] += [[items[1]] + [float(x) for x in items[2:]]]
            except KeyError:
                r[items[0]] = [[items[1]] + [float(x) for x in items[2:]]]
    if f == 1:
        ids = list(set([r[key][i][0] for key in list(r.keys()) for i in range(0, len(r[key]))]))
        return r, ids, f
    return r, list(set(taxids)), f


def process_lineages(items):
    """Looks up lineage information from taxids"""
    taxid, ranks, taxdir, dbname, lineage = items
    # Read the taxonomy db
    ncbi_taxa = init_sqlite_taxdb(taxdir, dbname)
    lineage_ranks = ncbi_taxa.get_rank(lineage)
    x = pd.DataFrame(lineage_ranks, index=["rank"]).T
    x = x.loc[x["rank"].isin(ranks)].reset_index().T
    x.columns = x.loc["rank"]
    x.drop("rank", inplace=True)
    x.index = [taxid]
    x = propagate_lower(x, taxid, ranks)
    x = add_names(x, taxid, ncbi_taxa)
    return x


def make_lineage_df(taxids, taxdir, dbname, ranks, threads=1):
    """Adds lineage information to diamond results"""
    # Read the taxonomy db
    ncbi_taxa = init_sqlite_taxdb(taxdir, dbname)
    lineages = ncbi_taxa.get_lineage_translator(taxids)
    items = [[taxid, ranks, taxdir, dbname, lineages[taxid]] for taxid in list(lineages.keys())]
    with Pool(processes=threads) as pool:
        res = list(
            tqdm.tqdm(pool.imap(process_lineages, items), desc="Creating lineages", total=len(items), unit=" taxids",
                      ncols=100))
    return pd.concat(res, sort=False)


def process_queries(items):
    """Receives a query and its results and assigns taxonomy"""
    res_tax = {}
    res_taxids = {}
    min_rank_threshold = 0
    query, lineage_df, res, rank_thresholds, args, taxidmap, lengthmap = items
    if len(rank_thresholds) > 0:
        min_rank_threshold = min([x for x in rank_thresholds.values()])
    columns = ['sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
               'sstart', 'send', 'evalue', 'bitscore']
    if len(res[0]) == 13:
        columns += ['staxids', 'slen']
    # Create pandas dataframe for slice
    res_df = pd.DataFrame(res, columns=columns, index=[query] * len(res))
    # Add taxidmap if not present in results
    if "staxids" not in res_df.columns:
        res_df = pd.merge(res_df, taxidmap, left_on="sseqid", right_index=True, how="left")
    # Add protein lengths if present
    if "slen" not in res_df.columns and len(lengthmap) > 0:
        if args.querylenmap:
            res_df = pd.merge(res_df, pd.DataFrame(lengthmap).T, left_index=True, right_index=True, how="left")
        elif args.subjectlenmap:
            res_df = pd.merge(res_df, lengthmap, left_on="sseqid", right_index=True, how="left")
    # Calculate normalized %id for slice
    res_df = calculate_norm_id(res_df, nolen=args.nolen)
    # Calculate bit score threshold for slice
    thresholds = get_thresholds(res_df, top=args.top)
    # Set index
    res_df.index.name = "qseqid"
    # Merge with lineage df
    res_df = pd.merge(res_df, lineage_df, left_on="staxids", right_index=True, how="left")
    # Initialize dictionaries
    res_tax[query] = dict.fromkeys(args.ranks, "Unclassified")
    res_taxids[query] = dict.fromkeys(args.ranks, -1)
    # Handle queries that return pandas Series
    res_df = res_df.loc[res_df.bitscore >= thresholds[query]]
    res_df = series2df(res_df)
    lca = {}
    lca_taxids = {}
    # Parse with rank thresholds or by just filtering by bitscore
    if "rank" in args.mode:
        if len(res_df.loc[res_df.normid >= min_rank_threshold]) > 0:
            lca, lca_taxids = parse_with_rank_thresholds(res_df, args.ranks, rank_thresholds, args.mode,
                                                         args.vote_threshold)
    else:
        lca, lca_taxids = get_lca(res_df, args.ranks)
    # Update results with lca and lca_taxids
    res_tax[query].update(lca)
    res_taxids[query].update(lca_taxids)
    res_tax[query] = propagate_taxa(res_tax[query], args.ranks)
    return res_tax[query], res_taxids[query], query


def write_blobout(f, res_taxids, queries, ranks):
    """Writes blob-format output for use with blobtools"""
    rev_ranks = [ranks[x] for x in list(range(len(ranks) - 1, -1, -1))]
    with open(f, 'w') as fhout:
        for i, query in enumerate(queries):
            d = res_taxids[i]
            for rank in rev_ranks:
                if rank in d.keys():
                    taxid = d[rank]
                    if taxid != -1:
                        fhout.write("{query}\t{taxid}\t1\tref\n".format(query=query, taxid=taxid))
                        break


def parse_hits(args):
    """
    Main function to handle diamond result files.

    Results are first read into a dictionary with queries as keys and values being a nested list with each item being
    the

    :param args: Input arguments from __main__.py
    :return: Return code 0 if all goes well
    """
    # Set up rank thresholds
    if "rank" in args.mode:
        rank_thresholds = get_rank_thresholds(args)
    else:
        rank_thresholds = {}
    # Read diamond results
    res, ids, f = read_df(args.diamond_results, args.taxidmap)
    # Read protein -> taxidmap file if specified
    taxidmap = pd.DataFrame()
    lengths = pd.DataFrame()
    if f == 1:
        taxidmap, taxids = read_taxidmap(args.taxidmap, ids)
        if args.querylenmap or args.subjectlenmap:
            lengthmap = read_lengthmap(args, res)
    else:
        taxids = ids
    # Create lineage dataframe
    lineage_df = make_lineage_df(taxids, args.taxdir, args.sqlitedb, args.ranks, args.threads)
    # Set up multiprocessing pool
    total_queries = len(res)
    items = []
    for q in res.keys():
        # If the diamond output does not have standard tango format (i.e. contains subject length and taxid) we
        # do some work to add this information.
        if f == 1:
            # Get all subject ids
            s = list(set([res[q][i][0] for i in range(0, len(res[q]))]))
            # Get lengths if specified by user
            if args.querylenmap:
                lengths = lengthmap.loc[q]
            elif args.subjectlenmap:
                lengths = lengthmap.loc[s]
            items.append([q, lineage_df, res[q], rank_thresholds, args, taxidmap.loc[s], lengths])
        # If diamond output has both taxonomy id and length then directly create the list of results to
        # feed into the multiprocessing pool
        else:
            items.append([q, lineage_df, res[q], rank_thresholds, args, None, None])
    with Pool(processes=args.threads) as pool:
        res = list(tqdm.tqdm(pool.imap(process_queries, items, chunksize=args.chunksize), desc="Parsing queries",
                             total=total_queries, unit=" queries", ncols=100))
    # res_tax is the taxonomy table with taxon names
    res_tax = [item[0] for item in res]
    # res_taxids is the taxonomy table with taxon ids
    res_taxids = [item[1] for item in res]
    # queries is a list of queries
    queries = [item[2] for item in res]
    # Writes blobtools-compatible output
    if args.blobout:
        sys.stderr.write("Writing blobtools file to {}\n".format(args.blobout))
        write_blobout(args.blobout, res_taxids, queries, args.ranks)
    # Create dataframe from taxonomy table
    res_df = pd.DataFrame(res_tax, index=queries)[args.ranks]
    res_df.index.name = "query"
    # Write main output
    sys.stderr.write("Writing main output to {}\n".format(args.outfile))
    res_df.to_csv(args.outfile, sep="\t")
    # Summary stats
    unc = [len(res_df.loc[res_df.loc[:, rank].str.contains("Unclassified")]) for rank in args.ranks]
    tot = [len(res_df)] * len(args.ranks)
    cl = 100 - np.divide(unc, tot) * 100
    cl = ["{}%".format(str(np.round(x, 1))) for x in cl]
    summary = pd.DataFrame(cl, index=args.ranks)
    sys.stderr.write("### SUMMARY ###:\nClassified sequences per rank:\n")
    summary.to_csv(sys.stderr, sep="\t", header=False)
    sys.stderr.write("\n")
    return 0
