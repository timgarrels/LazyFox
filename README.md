# LazyFox: Fast and parallelized overlapping community detection in large graphs

LazyFox is a parallelization enabled implementation of the Fox algorithm presented by Lyu et al. [1]

It allows overlapping community detection in very large graph datasets with billions of edges.

## Requirements
The LazyFox release binary has no external requirements.

To run the python code and notebooks you will need to have the following packages:
- pandas
- matplotlib
- [networkit](https://github.com/networkit/networkit) (you need to have cmake, a c++ compiler and Cython installed)
- jupyter (to run the notebooks)

## Usage

For a detailed step by step usage guide, refer to the jupyter notebooks in `notebooks/`.

Quick Instructions:

1. Download the latest release of LazyFox or compile it yourself.
2. Download a dataset to run on. In our project we used 
[Eu-core](https://snap.stanford.edu/data/email-Eu-core.html),
[DBLP](https://snap.stanford.edu/data/com-DBLP.html),
[LiveJournal](https://snap.stanford.edu/data/com-LiveJournal.html) and the
[Friendster](https://snap.stanford.edu/data/com-Friendster.html) datasets provided by the SNAP group.

You can use the `download.py` utility script for downloading them:
```
usage: download.py [-h] --dataset {eu,dblp,lj} --output OUTPUT

Download and Extract SNAP Datasets

arguments:
  -h, --help            show this help message and exit
  --dataset {eu,dblp,lj}
                        The dataset to download
  --output OUTPUT       The directory to download to
```
(Note that we do not provide automated download of the Friendster dataset due to its size.)

3. Run LazyFox from the commandline using:
```
usage: ./LazyFox input_graph output_directory [queue_size] [thread_count] [dumping {0,1}]

Run the LazyFox algorithm

arguments:
  input_graph         File containing the dataset as an edge list
  output_directory    Directory to log into, a subdirectory will be created for each run

optional arguments:
  queue_size          Specify the degree of parallel processing
  thread_count        Specify how many threads to use. Should be below or equal to queue_size
  dumping             If set to 0, LazyFox will not save clusterToNode results to disk.
  clustering          If set to 0, LazyFOX will use an external clustering. If this option is used, provide the filename of the clustering behind.
                      Every line of the file is assumed to consist of cluster_id node_id node_id node_id ... and no node_id can be present in multiple initial clusters.
  post-processing     If set to 0, LazyFOX will use an external script for postprocessing. If this option is used, provide the script path and paramters behind, lazy fox will invoke the script by adding the resulting file path behind.
```


Be sure to have enough cores to support the specified degree of parallelization.

## Example
```bash
# Download the DBLP dataset to the 'datasets/' directory
$ python3 download.py --dataset dblp --output datasets
# Make the LazyFox binary executable
$ chmod +x LazyFox
# Run LazyFox on the dblp graph, save the results to the 'output/' directory using a queue size of 2 and a thread count of 2
$ ./LazyFox datasets/com-dblp.ungraph.txt output 2 2
# Run LazyFox on the dblp graph, save the results to the 'output/' directory using a queue size of 2 and a thread count of 2, enabled dumping and an (initiaal) clustering from the file cluster_dblp.txt and no external postprocessing
$ ./LazyFox datasets/email-Eu-core.txt output 2 2 1 0 data/cluster_dblp.txt
# Run LazyFox on the dblp graph, save the results to the 'output/' directory using a queue size of 2 and a thread count of 2, enabled dumping, the internal clustering algorithm and invoke the user provided evaluate.py with the result file name as parameter once the algorithm has run completely.
$ ./LazyFox datasets/email-Eu-core.txt output 2 2 1 1 0 evaluate.py
```

[1] https://www.researchgate.net/publication/343408917_FOX_Fast_Overlapping_Community_Detection_Algorithm_in_Big_Weighted_Networks
