import urllib.request
import gzip
from os.path import join, exists
import shutil

urls = {
    "eu": [
        "https://snap.stanford.edu/data/email-Eu-core.txt.gz",
        "https://snap.stanford.edu/data/email-Eu-core-department-labels.txt.gz",
    ],
    "dblp": [
        "https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz",
        "https://snap.stanford.edu/data/bigdata/communities/com-dblp.all.cmty.txt.gz",
    ],
    "lj": [
        "https://snap.stanford.edu/data/bigdata/communities/com-lj.ungraph.txt.gz",
        "https://snap.stanford.edu/data/bigdata/communities/com-lj.all.cmty.txt.gz",
    ],
}


def download(dataset, directory):
    """Downloads and extracts dataset to the specified directory"""
    if dataset not in urls.keys():
        raise ValueError(f"No urls known for dataset {dataset}")
    for url in urls[dataset]:
        gz_filename = url.split("/")[-1]
        gz_path = join(directory, gz_filename)
        if not exists(gz_path):
            urllib.request.urlretrieve(url, gz_path)

        txt_filename = gz_filename.replace(".gz", "")
        txt_path = join(directory, txt_filename)
        if not exists(txt_path):
            with gzip.open(gz_path, "rb") as f_in:
                with open(txt_path, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)