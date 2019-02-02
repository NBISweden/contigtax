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

from ete3 import NCBITaxa
import pandas as pd
import subprocess
import os
import sys


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
    thresholds = (df.sort_values("bitscore", ascending=False).groupby(level=0).first().bitscore*((100-top))/100).to_dict()
    return thresholds


def add_names(x, taxid, ncbi_taxa):
    names = ncbi_taxa.get_taxid_translator(list(x.loc[taxid].values)+[taxid])
    n = {}
    for rank in list(x.columns):
        t = x.loc[taxid,rank]
        if t<0:
            name = "{}.{}".format(names[-t],rank)
        else:
            name = names[t]
        n["{}.name".format(rank)] = name
    n["species.name"] = names[taxid]
    name_df = pd.DataFrame(n, index=[taxid])
    return pd.merge(x, name_df, left_index=True, right_index=True)


def propagate_lower(x,taxid,ranks):
    rev_ranks = [ranks[x] for x in list(range(len(ranks) - 1, -1, -1))]
    missing = {}
    known = taxid
    for rank in rev_ranks[1:]:
        if not rank in x.columns:
            missing[rank] = -known
        else:
            known = x.loc[taxid,rank]
    return pd.merge(x, pd.DataFrame(missing,index=[taxid]), left_index=True, right_index=True)


def get_lineage_df(res, ncbi_taxa, ranks):
    lineages = ncbi_taxa.get_lineage_translator(res.staxids.unique())
    sys.stderr.write("#Extracting lineages for {} taxids\n".format(len(lineages)))
    lineage_df = pd.DataFrame(columns=set(ranks).difference(["species"]))
    i = 0
    for taxid, lineage in lineages.items():
        i+=1
        if i%100 == 0:
            p = float(i) / len(lineages)*100
            sys.stderr.write("#PROGRESS:   {}% ({}/{})\r".format(round(p,2), i, len(lineages)))
            sys.stderr.flush()
        lineage_ranks = ncbi_taxa.get_rank(lineage)
        x = pd.DataFrame(lineage_ranks,index=["rank"]).T
        x = x.loc[x["rank"].isin(ranks)].reset_index().T
        x.columns = x.loc["rank"]
        x.drop("rank", inplace=True)
        x.set_index("species", inplace=True)
        if len(x.columns)<len(ranks)-1:
            x = propagate_lower(x, taxid, ranks)
        x = add_names(x, taxid, ncbi_taxa)
        lineage_df = pd.concat([lineage_df, x], sort=True)
    lineage_df.index.name = "species"
    lineage_df = lineage_df.assign(species=pd.Series(lineage_df.index, index=lineage_df.index))
    sys.stderr.write("#PROGRESS:   100% ({}/{})\r".format(i,i))
    sys.stderr.flush()
    sys.stderr.write("\n")
    lineage_df[ranks].to_csv("lineage_df.tsv", sep="\t")
    return lineage_df


def get_lca(r, ranks):
    """Finds the rank where there's only one unique taxid and returns higher taxids"""
    query = r.index.unique()[0]
    rev_ranks = [ranks[x] for x in list(range(len(ranks) - 1, -1, -1))]
    for i, rank in enumerate(rev_ranks):
        higher_ranks = rev_ranks[i:]
        higher_rank_names = ["{}.name".format(x) for x in higher_ranks]
        c = r.groupby(rank).count()
        if len(c) == 1:
            if len(r)==1:
                lca = r.loc[query, higher_rank_names].values
                lca_taxids = r.loc[query, higher_ranks].values
            else:
                lca = r.loc[query, higher_rank_names].values[0]
                lca_taxids = r.loc[query, higher_ranks].values[0]
            return dict(zip(higher_ranks, lca)), dict(zip(higher_ranks, lca_taxids))
    return {}


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
    rank_df = r.groupby("staxids").first().reset_index()[[rank,"staxids"]].set_index("staxids")
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
            if known != "":
                res[rank] = "{}.{}".format("Unclassified",known)
        else:
            known = res[rank]
    return res


def series2df(df):
    if str(type(df)) == "<class 'pandas.core.series.Series'>":
        df = pd.DataFrame(df).T
    return df


def parse_hits6(res, top, ranks, mode, rank_thresholds, vote_threshold):
    """Parse the typical outfmt 6 blast file"""
    thresholds = get_thresholds(res, top=top)
    res_tax = {}
    res_taxids = {}
    queries = res.index.unique()
    if "rank" in mode:
        # Calculate normalized id values for each hit if rank thresholds are given
        sys.stderr.write("#Calculating normalized %id for hits\n")
        res = calculate_norm_id(res)
        min_rank_threshold = min([x for x in rank_thresholds.values()])
    else:
        queries = res.index.unique()
    sys.stderr.write("#Parsing blast file\n")
    for i, query in enumerate(queries, start=1):  # type: (int, object)
        if i % 10 == 0:
            p = (float(i)/len(queries))*100
            sys.stderr.write("#PROGRESS:   {}% ({}/{})\r".format(round(p, 2), i, len(queries)))
            sys.stderr.flush()
        # Initialize results for query
        res_tax[query] = dict.fromkeys(ranks, "Unclassified")
        res_taxids[query] = dict.fromkeys(ranks, -1)
        # Filter results for query
        _r = res.loc[query]
        _r = series2df(_r)
        _r = _r.loc[_r.bitscore >= thresholds[query]]
        _r = series2df(_r)
        if "rank" in mode:
            if len(_r.loc[_r.normid >= min_rank_threshold]) == 0:
                continue
            lca, lca_taxids = parse_with_rank_thresholds(_r, ranks, rank_thresholds, mode, vote_threshold)
        else:
            lca, lca_taxids = get_lca(_r, ranks)
        res_tax[query].update(lca)
        res_taxids[query].update(lca_taxids)
        res_tax[query] = propagate_taxa(res_tax[query], ranks)
    sys.stderr.write("#PROGRESS:   100% ({}/{})\n".format(i,len(queries)))
    return res_tax, res_taxids


def parse_hits102(taxids, ncbi_taxa,
                  ranks=["superkingdom","kingdom","phylum","class","order","family","genus","species"]):
    '''Parses the -102 format output from diamond'''
    lca_dict = {}
    for t in taxids:
        lineage = ncbi_taxa.get_lineage(t)
        lca_dict[t] = {}
        for rank in ranks:
            lca_dict[t][rank] = "Unclassified"
        for taxid in lineage:
            rank = ncbi_taxa.get_rank([taxid])[taxid]
            name = ncbi_taxa.get_taxid_translator([taxid])[taxid]
            if rank in ranks:
                lca_dict[t][rank] = name
        # Propagate lowest known taxonomy
        lca_dict[t] = propagate_taxa(lca_dict[t], ranks)
    return lca_dict


def init_taxdb(dbfile):
    if not dbfile:
        ncbi_taxa = NCBITaxa()
    else:
        dirname = os.path.dirname(dbfile)
        cmd = "mkdir -p {}".format(dirname)
        p1 = subprocess.Popen(cmd, shell=True, stdin=None)
        p1.wait()
        cmd = "touch {}".format(dbfile)
        p2 = subprocess.Popen(cmd, shell=True, stdin=None)
        p2.wait()
        ncbi_taxa = NCBITaxa(dbfile)
    return ncbi_taxa


def write_blobout(f, res_taxids, ranks):
    rev_ranks = [ranks[x] for x in list(range(len(ranks) - 1, -1, -1))]
    with open(f, 'w') as fhout:
        for query, d in res_taxids.items():
            for rank in rev_ranks:
                if rank in d.keys():
                    taxid = d[rank]
                    fhout.write("{query}\t{taxid}\t1\tref\n".format(query=query, taxid=taxid))
                    break