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
from multiprocessing import Pool
import tqdm
import sys
import gzip as gz
from tango.prepare import init_sqlite_taxdb


def calculate_norm_id(df):
    """Calculates normalized percent id for each hit"""
    # Create a new column that is: <alignment length>/<subject length> * <alignment identity>/100
    df = df.assign(lfrac=df.length.div(df.slen))
    # For hits where alignment length is greater than subject length, take the inverse
    df.loc[df.lfrac > 1, "lfrac"] = 1 / df.loc[df.lfrac > 1, "lfrac"]
    df = df.assign(normid=df.lfrac.multiply(df.pident.div(100)))
    df.drop("lfrac", axis=1, inplace=True)
    return df


def get_thresholds(df, top=10):
    """Get score threshold for each query"""
    thresholds = (df.sort_values("bitscore", ascending=False).groupby(level=0).first().bitscore * (
        (100 - top)) / 100).to_dict()
    return thresholds


def get_rank_thresholds(args):
    """Checks rank thresholds and creates dictionary to hold them"""
    rank_thresholds = {}
    assert (len(args.rank_thresholds) == len(args.ranks))
    for i, rank in enumerate(args.ranks):
        rank_thresholds[rank] = args.rank_thresholds[i]
    return rank_thresholds


def add_names(x, taxid, ncbi_taxa):
    """Adds taxa names for taxonomy ids"""
    names = ncbi_taxa.get_taxid_translator(list(x.loc[taxid].values)+[taxid])
    n = {}
    for rank in list(x.columns):
        t = x.loc[taxid,rank]
        if t<0:
            known_name = names[-t]
            if known_name == "root":
                name = "Unclassified"
            else:
                name = "{}.{}".format("Unclassified",known_name)
        else:
            name = names[t]
        n["{}.name".format(rank)] = name
    name_df = pd.DataFrame(n, index=[taxid])
    return pd.merge(x, name_df, left_index=True, right_index=True)


def propagate_lower(x, taxid, ranks):
    """Propagates taxonomy to lower ranks"""
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
    """Finds the rank where there's only one unique taxid and returns higher taxids"""
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


def read_df(f, open_function):
    """Reads the blast output into a format that can be processed with the multiprocessing module"""
    r = {}
    taxids = []
    with open_function(f, 'rt') as fhin:
        for i, line in enumerate(tqdm.tqdm(fhin, desc="Reading {}".format(f), ncols=100, unit=" lines")):
            line = line.rstrip()
            items = line.rsplit()
            taxids.append(items[-2])
            try:
                r[items[0]] += [[items[1]] + [float(x) for x in items[2:]]]
            except KeyError:
                r[items[0]] = [[items[1]] + [float(x) for x in items[2:]]]
    return r, list(set(taxids))


def process_lineages(items):
    """Looks up lineage information from taxids"""
    taxid, ranks, taxdir, lineage = items
    # Read the taxonomy db
    ncbi_taxa = init_sqlite_taxdb(taxdir)
    lineage_ranks = ncbi_taxa.get_rank(lineage)
    x = pd.DataFrame(lineage_ranks, index=["rank"]).T
    x = x.loc[x["rank"].isin(ranks)].reset_index().T
    x.columns = x.loc["rank"]
    x.drop("rank", inplace=True)
    x.index = [taxid]
    x = propagate_lower(x, taxid, ranks)
    x = add_names(x, taxid, ncbi_taxa)
    return x


def make_lineage_df(taxids, taxdir, ranks, threads=1):
    """Adds lineage information to diamond results"""
    # Read the taxonomy db
    ncbi_taxa = init_sqlite_taxdb(taxdir)
    lineages = ncbi_taxa.get_lineage_translator(taxids)
    items = [[taxid, ranks, taxdir, lineages[taxid]] for taxid in list(lineages.keys())]
    with Pool(processes=threads) as pool:
        res = list(tqdm.tqdm(pool.imap(process_lineages, items), desc="Creating lineages", total=len(items), unit=" taxids", ncols=100))
    return pd.concat(res, sort=False)


def process_queries(items):
    """Receives a query and its results and assigns taxonomy"""
    res_tax = {}
    res_taxids = {}
    min_rank_threshold = 0
    query, lineage_df, res, rank_thresholds, args = items
    if len(rank_thresholds) > 0:
        min_rank_threshold = min([x for x in rank_thresholds.values()])
    # Create pandas dataframe for slice
    res_df = pd.DataFrame(res, columns=['sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                                        'sstart', 'send', 'evalue', 'bitscore', 'staxids', 'slen'],
                          index=[query] * len(res))
    # Calculate normalized %id for slice
    res_df = calculate_norm_id(res_df)
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
            lca, lca_taxids = parse_with_rank_thresholds(res_df, args.ranks, rank_thresholds, args.mode, args.vote_threshold)
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
    """Parse the diamond blastx output"""
    # Set up rank thresholds
    if "rank" in args.mode:
        rank_thresholds = get_rank_thresholds(args)
    else:
        rank_thresholds = {}
    # Read diamond results
    if (args.diamond_results).split(".")[-1] == "gz":
        open_function = gz.open
    else:
        open_function = open
    res, taxids = read_df(args.diamond_results, open_function)
    # Create lineage dataframe
    lineage_df = make_lineage_df(taxids, args.taxdir, args.ranks, args.threads)
    # Set index to object
    #lineage_df.rename(index=lambda x: str(x), inplace=True)
    # Set up multiprocessing pool
    total_queries = len(res)
    items = [[q, lineage_df, res[q], rank_thresholds, args] for q in res.keys()]
    with Pool(processes=args.threads) as pool:
        res = list(tqdm.tqdm(pool.imap(process_queries, items, chunksize=args.chunksize), desc="Parsing queries", total=total_queries, unit=" queries", ncols=100))
    res_tax = [item[0] for item in res]
    res_taxids = [item[1] for item in res]
    queries = [item[2] for item in res]
    # Writes blobtools-compatible output
    if args.blobout:
        sys.stderr.write("Writing blobtools file to {}\n".format(args.blobout))
        write_blobout(args.blobout, res_taxids, queries, args.ranks)
    res_df = pd.DataFrame(res_tax, index=queries)[args.ranks]
    res_df.index.name = "query"
    sys.stderr.write("Writing main output to {}\n".format(args.outfile))
    res_df.to_csv(args.outfile, sep="\t")
    return 0