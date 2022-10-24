# LazyFox: Fast and parallelized overlapping community detection in large graphs

[LazyFox](https://arxiv.org/abs/2210.03211)[1] is a parallelization enabled implementation of the 
[FOX algorithm](https://dl.acm.org/doi/10.1145/3404970) presented by Lyu et al[2].

It allows overlapping community detection in very large graph datasets with billions of edges by optimizing a WCC estimation.
## Requirements
The LazyFox release binary has no external requirements. LazyFox can be compiled from source via cmake.

To run the python code and notebooks require the following packages:
- pandas
- matplotlib
- networkit (cmake, a c++ compiler and Cython required)
- jupyter (to run the notebooks)

## Quick Start
```bash
# Download the DBLP dataset to the 'datasets/' directory
$ python3 download.py --dataset dblp --output datasets
# Make the LazyFox binary executable
$ chmod +x LazyFox
# Run LazyFox on the dblp graph, save the results to the 'output/' directory using a queue size of 2 and a thread count of 2
$ ./LazyFox --input-graph datasets/com-dblp.ungraph.txt --output-dir output --queue-size 2 --thread-count 2
```

## Overview
The following section describes which input formats LazyFox accepts and how to interpret the output.

### Usage
For a detailed step-by-step usage guide, refer to the jupyter notebooks in `notebooks/`.
```
Run the LazyFox algorithm

Overlapping community detection for very large graph datasets with billions of edgesby optimizing a WCC estimation.

usage: LazyFox --input-graph <dataset path> --output-dir <output directory>
usage: LazyFox --input-graph <dataset path> --output-dir <output directory> --queue-size 64 --thread-count 64
usage: LazyFox --input-graph <dataset path> --output-dir <output directory> --wcc-threshold 0.05
usage: LazyFox --input-graph <dataset path> --output-dir <output directory> --pre-clustering <clustering path> --post-processing <script path>

Required Arguments:
  --input-graph <file_path>         File containing the graph dataset as an edge list
  --output-dir <directory>          Output directory, a subdirectory will be created for each run

Optional Arguments:
  --queue-size <int>                The degree of parallel processing (default=1)
  --thread-count <int>              How many threads to use. Should be below or equal to queue_size (default=1)
  --wcc-threshold <float>           Threshold in wcc-change to stop processing (default=0.01)
  --disable-dumping                 LazyFox will only save the clustering results of the final iteration to disk, not intermediate results
  --pre-clustering <file_path>      Loads external node clustering, replacing initial clustering algorithm
  --post-processing <script_path>   Script will be called after LazyFox clustering, with the run-subdirectory as argument
```

### Graph Input
LazyFox takes an edge list representing the initial graph as input file:
```
# Edge List for the following graph:
#    1     4     7
#  /   \ /   \ /   \
# 0     3     6     9
#  \   / \   / \   /
#    2     5     8
# 
# SourceNode	TargetNode
0   1
0   2
1   3
2   3
3   4
3   5
4   6
5   6
6   7
6   8
7   9
8   9
```
Each line contains one edge, consisting of the source node and target node, separated by a space or a tab.
The graph is undirected, therefore a single declaration of an edge is sufficient
(if node 0 and 1 are connected, only the line `0 1` is necessary, `1 0` is not).

### Pre-Clustering
LazyFox can work starting from a pre-existing node clustering. A pre-clustering can be supplied with the `--pre-clustering` argument: 
```
# Pre-Clustering for the nodes 3, 4, 5, and 6 of the graph above
3 4 5 6
```
The supplied file should contain one community per line, listing all nodes of that community in that line, separated by spaces or tabs.
Each node not mentioned in the pre-clustering is assumed to be a member of a single node community.

### Output

The output directory of LazyFox for the above example looks like this:
```
[output_directory]
    > CPP_[dataset_name]_[Timestamp]
        > iterations
            > 0.json
            > 0clusters.txt
```
The number in the filename states the iteration. `0.json` describes meta-information about the result of iteration 0, while `0clusters.txt` lists the generated communities in that run.
```
# 0.json
{
    "change_counter": {
        "copy": 0,
        "leave": 0,
        "stay": 10,
        "transfer": 0
    },
    "iteration": 0,
    "runtime": 7.5735e-05,
    "wcc_lookup": {
        "0": null
    }
}
```
```
# 0clusters.txt
3	4	5	6
```
The `json` file lists the iterations runtime, wcc-scores per community, and the node changes performed during the iteration.
The `txt` file lists one community per line, each line containing the node ids for that community.

Due to the small size of this example, the original community was preserved by LazyFox. Larger graphs will result in multiple iterations and multiple `json` and `txt` files.

### Download Utility
You can use the `download.py` utility script to download the SNAP group datasets we used:
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

## References
If you use LazyFox, please cite
 
[1] Garrels, T., Khodabakhsh, A., Renard, B. Y. & Baum, K. LazyFox: Fast and parallelized overlapping community detection in large graphs. doi:10.48550/ARXIV.2210.03211 (2022)

[2] Lyu, T., Bing, L., Zhang, Z. & Zhang, Y. FOX: Fast Overlapping Community Detection Algorithm in Big Weighted Networks. Trans. Soc. Comput. 3, Article 16, doi:10.1145/3404970 (2020)
