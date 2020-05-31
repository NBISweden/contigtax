import pytest
import pandas as pd
from tango.assign import (get_rank_thresholds, get_lca, get_thresholds,
                          propagate_taxids, propagate_lower)


@pytest.mark.parametrize("top,threshold", [(10, 90), (5, 95)])
def test_get_thresholds(top, threshold):
    df = pd.DataFrame({"sseqid": ["subject{}".format(i) for i in range(1, 6)],
                       "pident": [90, 85, 80, 75, 70],
                       "evalue": [0.0001, 0.001, 0.01, 0.1, 1],
                       "bitscore": [100, 90, 80, 70, 60]}, index=["query"] * 5)
    assert get_thresholds(df, top)["query"] == threshold


def test_get_rank_thresholds():
    thresholds = get_rank_thresholds(ranks=["phylum", "genus", "species"],
                                     thresholds=[45, 60, 80])
    assert thresholds["phylum"] == 45
    assert thresholds["genus"] == 60
    assert thresholds["species"] == 80


def test_lca():
    assignranks = ["phylum", "genus", "species"]
    reportranks = ["phylum", "genus", "species"]
    r = pd.DataFrame(
        {"species": [1, 2, 3], "genus": [11, 11, 22], "phylum": [111, 111, 111],
         "superkingdom": [2, 2, 2]}, index=["q1"] * 3)
    lca = get_lca(r, assignranks, reportranks)
    assert lca["phylum"] == 111


@pytest.mark.parametrize("res,taxid", [({'species': -1, 'family': -171549,
                                         'genus': -171549, 'order': 171549,
                                         'phylum': 976, 'class': 200643,
                                         'superkingdom': 2}, -171549), (
                                        {'species': -1, 'family': -2,
                                         'genus': -2, 'order': -2,
                                         'phylum': -2, 'class': -2,
                                         'superkingdom': 2}, -2)])
def test_propagate_taxids(res, taxid):
    ranks = ["superkingdom", "phylum", "class", "order", "family", "genus",
             "species"]
    assert propagate_taxids(res, ranks)["species"] == taxid


@pytest.mark.parametrize("lower,propagated", [(2, -2), (1117, -1117)])
def test_propagate_lower(lower, propagated):
    x = pd.DataFrame({lower: {"rank": lower}}).T
    ranks = ["rank", "superkingdom", "phylum", "class", "order", "family", "genus",
             "species"]
    assert propagate_lower(x, lower, ranks).loc[lower, "species"] == propagated