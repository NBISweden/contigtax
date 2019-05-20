import pandas as pd
import numpy as np
from multiprocessing import Pool
import tqdm
import sys
import gzip as gz
from tango.prepare import init_sqlite_taxdb


def get_thresholds(df, top=10):
    """
    Here bit-score thresholds are calculated per query an returned in a dictionary.

    The pandas DataFrame is first sorted by bitscore (high to low), then grouped by query, then for the first entry
    per query the top% of the best hit is calculated and converted to dictionary.

    Parameters
    ----------
    df: pandas.DataFrame
        DataFrame of diamond results
    top: int
        Percentage range of top bitscore

    Returns
    -------
    thresholds: dict
        Dictionary with queries as keys and bitscore thresholds as values
    """

    thresholds = (df.sort_values("bitscore", ascending=False).groupby(level=0).first().bitscore * (
        (100 - top)) / 100).to_dict()
    return thresholds


def get_rank_thresholds(ranks, thresholds):
    """
    Constructs dictionary of rank-specific thresholds

    Parameters
    ----------
    ranks: list
        Taxonomic ranks to assign
    thresholds: list
        Thresholds for taxonomic ranks

    Returns
    -------
        Dictionary of thresholds
    """

    t_len, r_len = len(thresholds), len(ranks)
    if t_len != r_len:
        sys.exit("ERROR: Number of taxonomic ranks ({}) and number of thresholds ({}) differ\n".format(r_len, t_len))
    return dict(zip(ranks, thresholds))


def add_names(x, taxid, ncbi_taxa):
    """
    This function translates taxonomy ids to names. It operates per-row in the lineage dataframe.

    Parameters
    ----------
    x: pandas.DataFrame
        DataFrame of one taxid and its taxonomic ranks
    taxid: int
        Taxid being evaluated
    ncbi_taxa: ete3.ncbi_taxonomy.ncbiquery.NCBITaxa
        The ete3 sqlite database connection

    Returns
    -------
        The original DataFrame merged with the taxa names
    """

    # Get a names dictionary for all taxids in the row
    names = ncbi_taxa.get_taxid_translator(list(x.loc[taxid].values) + [taxid])
    n = {}
    # Iterate ranks
    for rank in list(x.columns):
        # Get taxid for the current rank
        t = x.loc[taxid, rank]
        # If taxid is negative it means that there is no classified taxonomy at this rank
        # Instead we get the last known name in the hierarchy. We can then use the negative values to translate into
        # the name with the "Unclassified." prefix.
        # If the name is 'root' we just use 'Unclassified'
        if t < 0:
            known_name = names[-t]
            if known_name == "root":
                name = "Unclassified"
            else:
                name = known_name
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

    Parameters
    ----------
    x: pandas.DataFrame
        DataFrame of one taxid and its taxonomic ranks
    taxid:  int
        Taxid being evaluated
    ranks: list
        Ranks used for assigning

    Returns
    -------
        pandas.DataFrame updated with missing ranks

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


def get_lca(r, assignranks, reportranks):
    """
    Assign lowest common ancestor from a set of taxids.

    Parameters
    ----------
    r: pandas.DataFrame
        Results for a single query, extracted from the main diamond results file
    assignranks: list
        Taxonomic ranks to assign taxonomy for
    reportranks: list
        Taxonomic ranks to report taxonomy for

    Returns
    -------
        a tuple of dictionaries with ranks as keys and taxa names/ids as values

    This function takes a query-slice of the diamond results after filtering by score (and rank-threshold if tango mode
    is 'rank_lca' or 'rank_vote'). It then iterates through each rank in reverse order checks how many unique taxids are
    found at that rank. If there's only one taxid
    """

    query = r.index.unique()[0]
    # Reverse ranks for iterating
    rev_ranks = [assignranks[x] for x in list(range(len(assignranks) - 1, -1, -1))]
    # Iterate through the assignranks
    for rank in rev_ranks:
        higher_ranks = reportranks[0:reportranks.index(rank) + 1]
        higher_rank_names = ["{}.name".format(x) for x in higher_ranks]
        # Count number of taxa at rank
        c = r.groupby(rank).count()
        # If there's only one taxa then we have found the LCA
        if len(c) == 1:
            if len(r) == 1:
                lca_taxids = r.loc[query, higher_ranks].values
            else:
                lca_taxids = r.loc[query, higher_ranks].values[0]
            return dict(zip(higher_ranks, lca_taxids))
    return {}


def parse_with_rank_thresholds(r, assignranks, reportranks, rank_thresholds, mode, vote_threshold):
    """Assigns taxonomy using rank_specific thresholds

    The ranks used to assign taxonomy are iterated in reverse (e.g. species, genus, phylum),
    at each rank results are filtered by the corresponding rank threshold,
    if no hits remain after filtering the next rank is evaluated,

    Then, if mode=='rank_lca', for remaining hits, a lowest common ancestor is calculated from all remaining taxids.

    However, if mode=='rank_vote', taxids are counted among the remaining hits and all results matching taxids
    that occur more than vote_threshold are used to determine the lowest common ancestor.

    If a taxonomy can be assigned at a rank, it is returned directly. If no taxonomy can be assigned at any of the
    ranks, empty results are returned.

    Parameters
    ----------
    r: pandas.DataFrame
        Dataframe slice for a query
    assignranks: list
        Taxonomic ranks used to assign taxonomy
    reportranks: list
        Taxonomic ranks at which taxonomy is reported
    rank_thresholds: dict
        Dictionary of rank_specific thresholds
    mode: str
        'rank_lca' or 'rank_vote'
    vote_threshold: float
        Cutoff used to filter out common taxids

    Returns
    -------
    tuple
        Dictionaries with taxonomy names and taxonomy ids at each rank
    """

    # Start from lowest rank
    rev_ranks = [assignranks[x] for x in list(range(len(assignranks) - 1, -1, -1))]
    for rank in rev_ranks:
        # Make sure that LCA is not set below current rank
        allowed_ranks = assignranks[0:assignranks.index(rank) + 1]
        # Get rank threshold
        threshold = rank_thresholds[rank]
        # Filter results by rank threshold
        try:
            _r = r.loc[r.pident >= threshold]
        except KeyError:
            continue
        if len(_r) == 0:
            continue
        lca_taxids = {}
        # After filtering, either calculate lca from all filtered taxids
        if mode == "rank_lca":
            lca_taxids = get_lca(_r, allowed_ranks, reportranks)
        # Or at each rank, get most common taxid
        elif mode == "rank_vote":
            vote = get_rank_vote(_r, rank, vote_threshold)
            if len(vote) > 0:
                lca_taxids = get_lca(vote, allowed_ranks, reportranks)
        if len(lca_taxids.keys()) > 0:
            return lca_taxids
    return {}


def get_rank_vote(r, rank, vote_threshold=0.5):
    """
    Filter results based on fraction of taxa

    Parameters
    ----------
    r: pandas.DataFrame
        Results for a single query, after filtering with bitscore and rank-specific thresholds
    rank: str
        Current rank being investigated
    vote_threshold: float
        Required fraction of hits from a single taxa in order to keep taxa

    Returns
    -------
        Filtered dataframe only containing taxa that meet vote_threshold

    Here taxa are counted among all hits remaining for a query after filtering using bitscore and rank-specific
    thresholds. Taxa are counted at a certain rank and counts are normalized. Hits belonging to taxa above
    vote_threshold are kept while others are filtered out.
    """

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


def propagate_taxids(res, ranks):
    """
    Transfer taxonomy ids to unassigned ranks based on best known taxonomy
    Example:

    {'species': -1, 'family': -171549, 'genus': -171549, 'order': 171549, 'phylum': 976, 'class': 200643, 'superkingdom': 2}

    should become
    {'species': -171549, 'family': -171549, 'genus': -171549, 'order': 171549, 'phylum': 976, 'class': 200643, 'superkingdom': 2}

    Parameters
    ----------
    res: dict
        Dictionary of ranks and taxonomy ids
    ranks: list
        Ranks to assign taxonomy to

    Returns
    -------
    res: dict
        Dictionary with updated taxonomy ids
    """

    known = -1
    for rank in ranks:
        # If not -1 (Unclassified) at rank, store assignment as known
        if res[rank] != -1:
            known = res[rank]
            continue
        # If -1 at rank (Unclassified), add the taxid with the '-' prefix
        if res[rank] == -1:
            res[rank] = -abs(known)
    return res


def series2df(df):
    """Converts pandas series to pandas dataframe"""
    if str(type(df)) == "<class 'pandas.core.series.Series'>":
        df = pd.DataFrame(df).T
    return df


def read_taxidmap(f, ids):
    """
    Reads the protein to taxid map file and stores mappings

    Parameters
    ----------
    f: str
        Input file with protein_id->taxid map
    ids: list
        Protein ids to store taxids for

    Returns
    -------
        Dictionary of protein ids to taxid and all unique taxids
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


def read_df(infile, top=10, e=0.001, input_format="tango", taxidmap=None):
    """
    Reads the blast results from file and returns a dictionary with query->results.

    Note that the input is assumed to be sorted by bitscore for each query. The first entry for a query is used to set
    the score threshold for storing hits for that query. So if a query has a bitscore of 100 and --top 10 is specified
    then we only store subsequent hits that have a bitscore of at least (100-0.1*100) = 90.

    Tango-formatted output contains two additional compared to the standard blast format 6:
    query1     subject1   93.6    47      3       0       146     6       79      125     8.5e-16 91.3    314295
    query1     subject2  100.0   44      0       0       137     6       484     527     2.5e-15 89.7    9347
    query2     subject3      53.5    241     84      2       645     7       15      255     1.3e-53 216.9   864142

    where the last column is the taxid of the subject.

    Otherwise the output may have the typical blast format 6 output.

    Parameters
    ----------
    infile: str
        Arguments from argument parser
    top: int
        Keep results within top% of best bitscore
    e: float
        Maximum allowed e-value to keep a hit.
    input_format: str
        Blast format. 'tango' if taxid for each subject is present in blast results, otherwise 'blast'
    taxidmap: str
        File mapping each subject id to a taxid

    Returns
    -------
    tuple
        The function returns a tuple with dictionary of query->results and
        unique taxonomy ids (if tango format) or unique subject ids
    """

    open_function = open
    if ".gz" in infile:
        open_function = gz.open
    r = {}
    taxids = []
    queries = {}
    with open_function(infile, 'rt') as fhin:
        for line in tqdm.tqdm(fhin, desc="Reading {}".format(infile), ncols=100, unit=" lines"):
            items = line.rstrip().rsplit()
            query, subject, pident, evalue, score = items[0], items[1], float(items[2]), \
                                                    float(items[10]), float(items[11])
            try:
                min_score = queries[query]['min_score']
            except KeyError:
                min_score = score * ((100 - top) / 100)
                queries[query] = {'min_score': min_score}
            if score < min_score or evalue > e:
                continue
            if input_format == "tango" and len(items) > 12:
                taxid = items[12]
                taxids.append(taxid)
                # TODO: Is there a way to skip storing the same taxid from a worse hit for the same query
            elif input_format == "blast" and len(items) == 12:
                taxid = ""
                if not taxidmap:
                    sys.exit(
                        "ERROR: Standard blast input detected with no protein -> taxid file specified (--taxidmap).")
            else:
                continue
            # Add results for query to dictionary
            try:
                r[query] += [[subject, pident, evalue, score, int(taxid)]]
            except KeyError:
                r[query] = [[subject, pident, evalue, score, int(taxid)]]
    # If this is blast format then we return all subject ids found
    if input_format == "blast":
        ids = list(set([r[key][i][0] for key in list(r.keys()) for i in range(0, len(r[key]))]))
        return r, ids
    # If this is tango format then return all taxids found
    return r, list(set(taxids))


def process_lineages(items):
    """
    Looks up lineage information from taxids.

    The lineage object is a list of taxonomic ids corresponding to the full lineage of a single taxid.
    """
    taxid, ranks, taxdir, dbname, lineage = items
    # Read the taxonomy db
    ncbi_taxa = init_sqlite_taxdb(taxdir, dbname)
    # Get ranks for each taxid in the lineage
    lineage_ranks = ncbi_taxa.get_rank(lineage)
    x = pd.DataFrame(lineage_ranks, index=["rank"]).T
    x = x.loc[x["rank"].isin(ranks)].reset_index().T
    x.columns = x.loc["rank"]
    x.drop("rank", inplace=True)
    x.index = [taxid]
    # Add taxids for lower ranks in the hierarchy
    x = propagate_lower(x, taxid, ranks)
    # Add names for taxids
    x = add_names(x, taxid, ncbi_taxa)
    return x


def make_name_dict(df, ranks):
    """
    Creates a dictionary of taxids to taxonomy names, including Unclassified ranks

    Parameters
    ----------
    df: pandas.DataFrame
        Lineage dataframe
    ranks: list
        Ranks to store names information for

    Returns
    -------
    name_dict: dict
        Name dictionary mapping taxonomy ids to names
    """

    name_dict = {}
    for rank in ranks:
        name_dict.update(dict(zip(df[rank].values, df["{}.name".format(rank)].values)))
        name_dict.update(dict(zip(-abs(df[rank]), "Unclassified." + df["{}.name".format(rank)])))
    name_dict[-1] = "Unclassified"
    return name_dict


def make_lineage_df(taxids, taxdir, dbname, ranks, cpus=1):
    """
    Creates a lineage dataframe with full taxonomic information for a list of taxids.

    Example:
    taxid   species phylum  genus   genus.name      phylum.name     species.name
    859655  305     1224    48736   Ralstonia       Proteobacteria  Ralstonia solanacearum
    387344  1580    1239    1578    Lactobacillus   Firmicutes      Lactobacillus brevis
    358681  1393    1239    55080   Brevibacillus   Firmicutes      Brevibacillus brevis

    Parameters
    ----------
    taxids: list
        List of taxonomic ids to obtain information for
    taxdir: str
        Path to directory holding taxonomic info
    dbname: str
        Name of ete3 sqlite database within taxdir
    ranks: list
        Ranks to store information for
    cpus: int
        Number of cpus to use

    Returns
    -------
    lineage_df: pandas.DataFrame
        Data Frame with full taxonomic info
    """
    # Read the taxonomy db
    ncbi_taxa = init_sqlite_taxdb(taxdir, dbname)
    lineages = ncbi_taxa.get_lineage_translator(taxids)
    # Store potential missing taxids and warn user
    missing_taxids = set([int(x) for x in taxids]).difference(lineages.keys())
    # Get possible translations for taxids that have been changed
    _, translate_dict = ncbi_taxa._translate_merged(list(set(taxids).difference(lineages.keys())))
    rename = {y: x for x, y in translate_dict.items()}
    # Update lineages with missing taxids
    lineages.update(ncbi_taxa.get_lineage_translator(translate_dict.values()))
    items = [[taxid, ranks, taxdir, dbname, lineages[taxid]] for taxid in list(lineages.keys())]
    with Pool(processes=cpus) as pool:
        res = list(
            tqdm.tqdm(pool.imap(process_lineages, items), desc="Making lineages", total=len(items),
                      unit=" taxids", ncols=100))
    lineage_df = pd.concat(res, sort=False)
    lineage_df.rename(index=rename, inplace=True)
    lineage_df.rename(index=lambda x: int(x), inplace=True)
    for rank in ranks:
        lineage_df[rank] = pd.to_numeric(lineage_df[rank])
    name_dict = make_name_dict(lineage_df, ranks)
    if len(missing_taxids) > 0:
        sys.stderr.write("#WARNING: Missing taxids found:\n")
        sys.stderr.write("#{}\n".format(",".join([str(x) for x in missing_taxids])))
        sys.stderr.write("#To fix this, you can try to update the taxonomy database using\n")
        sys.stderr.write("#tango download taxonomy --force\n")
    return lineage_df.loc[:,lineage_df.dtypes==int], name_dict


def process_queries(args):
    """Receives a query and its results and assigns taxonomy"""
    res_taxids = {}
    min_rank_threshold = 0
    query, res, rank_thresholds, top, reportranks, assignranks, mode, vote_threshold, lineage_df, taxidmap = args
    if len(rank_thresholds) > 0 and "rank" in mode:
        min_rank_threshold = min([x for x in rank_thresholds.values()])
    columns = ['sseqid', 'pident', 'evalue', 'bitscore']
    if len(res[0]) == 5:
        columns += ['staxids']
    # Create pandas dataframe for slice
    res_df = pd.DataFrame(res, columns=columns, index=[query] * len(res))
    # Add taxidmap if not present in results
    if "staxids" not in res_df.columns:
        res_df = pd.merge(res_df, taxidmap, left_on="sseqid", right_index=True, how="left")
    # Calculate bit score threshold for slice
    thresholds = get_thresholds(res_df, top=top)
    # Set index
    res_df.index.name = "qseqid"
    # Merge with lineage df
    res_df = pd.merge(res_df, lineage_df, left_on="staxids", right_index=True, how="left")
    # Remove potential nan rows created if the blast results have taxids that are missing from lineage_df
    res_df = res_df.loc[res_df[reportranks[0]] == res_df[reportranks[0]]]
    # Initialize dictionary
    res_taxids[query] = dict.fromkeys(reportranks, -1)
    # Handle queries that return pandas Series
    res_df = res_df.loc[res_df.bitscore >= thresholds[query]]
    res_df = series2df(res_df)
    lca_taxids = {}
    # Parse with rank thresholds or by just filtering by bitscore
    if "rank" in mode:
        if len(res_df.loc[res_df.pident >= min_rank_threshold]) > 0:
            lca_taxids = parse_with_rank_thresholds(res_df, assignranks, reportranks,
                                                         rank_thresholds, mode, vote_threshold)
    else:
        lca_taxids = get_lca(res_df, assignranks, reportranks)
    # Update results with lca_taxids
    res_taxids[query].update(lca_taxids)
    res_taxids[query] = propagate_taxids(res_taxids[query], reportranks)
    return res_taxids[query], query


def write_blobout(f, res_taxids, queries, ranks):
    """
    Writes output in a format for use with blobtools

    Parameters
    ----------
    f: str
        Outputfile
    res_taxids: list
        List of results for queries
    queries: list
        List of queries
    ranks: list
        Ranks to write results for
    """

    rev_ranks = [ranks[x] for x in list(range(len(ranks) - 1, -1, -1))]
    with open(f, 'w') as fhout:
        for i, query in enumerate(queries):
            d = res_taxids[i]
            for rank in rev_ranks:
                if rank in d.keys():
                    taxid = d[rank]
                    if taxid != -1:
                        fhout.write("{query}\t{taxid}\t1\tref\n".format(query=query, taxid=abs(taxid)))
                        break


def stage_queries(res, lineage_df, input_format="tango", rank_thresholds=[45, 60, 85], top=10, mode="rank_lca",
                  vote_threshold=0.5, assignranks=["phylum", "genus", "species"],
                  reportranks=["superkingdom", "phylum", "class", "order", "family", "genus", "species"],
                  taxidmap=None):
    """

    Parameters
    ----------
    res: dict
        Dictionary with queries as keys and a list of hits as values
    lineage_df: pandas.DataFrame
        Data frame of taxids and taxonomic information
    input_format: str
        'tango' or 'blast'
    rank_thresholds: list
        List of thresholds for ranks
    top: int
        Only evaluate results within <top> percent bitscore of best scoring hit
    mode: str
        'rank_lca' or 'rank_vote' for rank thresholds usage or 'score' to just filter by bitscore
    vote_threshold: float
        Cutoff used to filter out common taxids
    assignranks: list
        Ranks used to assign taxonomy
    reportranks: list
        Ranks to report taxonomy for (inferred from assignranks)
    taxidmap: dict
        Dictionary with subject ids as keys and taxids as values

    Returns
    -------
    items: list
        List of items to send to multiprocessing
    """

    items = []
    total_queries = len(res)
    for q in tqdm.tqdm(sorted(res.keys()), total=total_queries, unit=" queries", ncols=100, desc="Staging queries"):
        # If the diamond output does not have standard tango format we do some work to add this information.
        item = [q, res[q], rank_thresholds, top, reportranks, assignranks, mode, vote_threshold]
        if input_format == "blast":
            # Get all subject ids
            s = list(set([res[q][i][0] for i in range(0, len(res[q]))]).intersection(lineage_df.index))
            # Get all taxids for this query
            q_taxids = taxidmap.loc[s, "staxids"].unique()
            item += [lineage_df.loc[q_taxids], taxidmap.loc[s]]
        # If diamond output has taxonomy id then directly create the list of results to
        # feed into the multiprocessing pool
        else:
            # Get all taxids for query
            q_taxids = list(set([res[q][i][-1] for i in range(0, len(res[q]))]).intersection(lineage_df.index))
            item += [lineage_df.loc[q_taxids], None]
        items.append(item)
    return items


def parse_hits(diamond_results, outfile, taxidout=False, blobout=False, top=10, evalue=0.001, input_format="tango",
               taxidmap=False, mode="rank_lca", vote_threshold=0.5, assignranks=["phylum", "genus", "species"],
               reportranks=["superkingdom", "phylum", "class", "order", "family", "genus", "species"],
               rank_thresholds=[45, 60, 85], taxdir="./taxonomy/", sqlitedb="taxonomy.sqlite", chunksize=1, cpus=1):
    """
    This is the main function to handle diamond result files and assign taxonomy to queries.

    The function performs the following steps:
    1. Checks rank-specific thresholds
    2. Reads the diamond results file
    3. If required, maps subject ids to taxonomy ids
    4. Creates a dataframe of all unique taxonomy ids found for subjects and their taxa names for each rank
    5. Stages queries for multiprocessing
    6. Processes each query and returns it with assigned taxonomy
    7. Writes output to file

    Parameters
    ----------
    diamond_results: str
        Diamond results file
    outfile: str
        File to write results to
    taxidout: str
        If True, write results with taxonomic ids instead of names to file
    blobout: str
        If True, write output in blobtools format
    top: int
        Evaluate hits within this bitscore percent range of the best scoring hit
    evalue: float
        Filter hits with evalue larger than this
    input_format: str
        'tango' or 'blast' depending on whether the diamond results has subject taxids in the last column or not
    taxidmap: str
        Path to a file mapping subject ids to taxids (needed if input_format != 'tango')
    mode: str
        How to assign taxonomy: 'rank_lca' and 'rank_vote' use rank specific thresholds,
        'score' only filters by bitscore
    vote_threshold: float
        When using 'rank_vote' to assign taxonomy, this is the fraction of hits that must have the same taxid to
        assign a taxonomy at a rank
    assignranks: list
        Ranks used to assign taxonomy
    reportranks: list
        Ranks to report taxonomy for
    rank_thresholds: list
        Percent identity thresholds for assigning taxonomy
    taxdir: str
        Path to directory holding taxonomic information files
    sqlitedb: str
        Name of ete3 sqlite database within taxdir
    chunksize: int
        The size of chunks for the iterable being submitted to the process pool
    cpus: int
        The number of worker processes to use
    args:
        Input arguments from __main__.py

    Returns
    -------
        Return code 0 if function finished without issue
    """

    # Set up rank thresholds
    if "rank" in mode:
        rank_thresholds = get_rank_thresholds(assignranks, rank_thresholds)
    # Read diamond results
    res, ids = read_df(diamond_results, top, evalue, input_format, taxidmap)
    # Read protein -> taxidmap file if specified
    taxidmap = pd.DataFrame()
    if input_format == "blast":
        taxidmap, taxids = read_taxidmap(taxidmap, ids)
    else:
        taxids = ids
    # Create lineage dataframe
    lineage_df, name_dict = make_lineage_df(taxids, taxdir, sqlitedb, reportranks, cpus)
    # Set up multiprocessing pool
    items = stage_queries(res, lineage_df, input_format, rank_thresholds, top, mode, vote_threshold, assignranks,
                          reportranks, taxidmap)
    total_queries = len(res)
    with Pool(processes=cpus) as pool:
        assign_res = list(tqdm.tqdm(pool.imap(process_queries, items, chunksize=chunksize), desc="Parsing queries",
                                    total=total_queries, unit=" queries", ncols=100))
    # res_tax is the taxonomy table with taxids
    res_tax = [item[0] for item in assign_res]
    # queries is a list of queries
    queries = [item[1] for item in assign_res]
    # Create dataframe from taxonomy results
    res_df = pd.DataFrame(res_tax, index=queries)[reportranks]
    res_df.index.name = "query"
    # Writes blobtools-compatible output
    if blobout:
        sys.stderr.write("Writing blobtools file to {}\n".format(blobout))
        write_blobout(blobout, res_df, queries, reportranks)
    # Write table with taxonomy ids instead of taxon names
    if taxidout:
        sys.stderr.write("Writing results with taxids to {}\n".format(taxidout))
        res_df.to_csv(taxidout, sep="\t")
    # Write main output
    sys.stderr.write("Translating taxids to names\n")
    res_names_df = pd.concat([res_df[rank].astype("category").cat.rename_categories(name_dict) for rank in reportranks],
                             axis=1)
    sys.stderr.write("Writing main output to {}\n".format(outfile))
    res_names_df.to_csv(outfile, sep="\t")
    # Summary stats
    unc = [len(res_names_df.loc[res_names_df.loc[:, rank].str.contains("Unclassified")]) for rank in reportranks]
    tot = [len(res_names_df)] * len(reportranks)
    cl = 100 - np.divide(unc, tot) * 100
    cl = ["{}%".format(str(np.round(x, 1))) for x in cl]
    summary = pd.DataFrame(cl, index=reportranks)
    sys.stderr.write("### SUMMARY ###:\nClassified sequences per rank:\n")
    summary.to_csv(sys.stderr, sep="\t", header=False)
    sys.stderr.write("\n")
    return 0