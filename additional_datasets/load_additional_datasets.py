"""Loads the network datasets as nx.Graph and provides them in a name->nx.Graph lookup dict"""

import networkx as nx
from zipfile import ZipFile
import tarfile
import io


def _get_football() -> nx.Graph:
    with ZipFile("datasets/football/football.zip") as zf:
        zf.extract("football.gml", "datasets/football/extracted")

    return nx.read_gml("datasets/football/extracted/football.gml")


def _get_karate_club() -> nx.Graph:
    with tarfile.open('datasets/karate_club/download.tsv.ucidata-zachary.tar.bz2', 'r:bz2') as f:
        edgelist_file = f.extractfile('ucidata-zachary/out.ucidata-zachary').read()

    return nx.read_edgelist(io.BytesIO(edgelist_file), comments="%")


def _get_celegans_metabolic() -> nx.Graph:
    return nx.read_adjlist("datasets/celegans_metabolic/celegans_metabolic.graph.bz2")


def _get_lesmis() -> nx.Graph:
    with ZipFile("datasets/miserables/lesmis.zip") as zf:
        zf.extract("lesmis.gml", "datasets/miserables/extracted")

    return nx.read_gml("datasets/miserables/extracted/lesmis.gml")


def _get_stelzl_human() -> nx.Graph:
    with tarfile.open('datasets/stelzl_human/download.tsv.maayan-Stelzl.tar.bz2', 'r:bz2') as f:
        edgelist_file = f.extractfile('maayan-Stelzl/out.maayan-Stelzl').read()

    return nx.read_edgelist(io.BytesIO(edgelist_file), comments="%")


def _get_fidgeys_human() -> nx.Graph:
    with tarfile.open('datasets/fidgeys_human/download.tsv.maayan-figeys.tar.bz2', 'r:bz2') as f:
        edgelist_file = f.extractfile('maayan-figeys/out.maayan-figeys').read()

    return nx.read_edgelist(io.BytesIO(edgelist_file), comments="%")


graphs = {
    "football": _get_football(),
    "celegans_metabolic": _get_celegans_metabolic(),
    "lesmis": _get_lesmis(),
    "karate_club": _get_karate_club(),
    "fidgeys_human": _get_fidgeys_human(),
    "stelzl_human": _get_stelzl_human(),
}
