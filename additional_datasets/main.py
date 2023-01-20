import os
from tqdm import tqdm
import networkx as nx
import json

from load_additional_datasets import graphs

LAZYFOX_BINARY = "../LazyFox"
OUTDIR = "output"
os.mkdir(OUTDIR)

# Technically, the maximum queue sizes should be scaled by the size of the respective dataset
# For example, sqrt(node_count), or more complex heuristics like queue_size * avg_node_degree >= X% * edge_count
# However, our largest dataset currently accounts for 2239 nodes, and 10% of that includes all such business
# for all datasets, including itself
QUEUE_SIZES = list(range(1, 301))


def run_lazyfox(
        binary_path: str,
        input_graph_path: str,
        output_directory: str,
        queue_size: int,
        thread_count: int,
        dumping: bool,
        log_file: str = "std.log",
        benchmark_file: str = "bench.mark",
        config_file: str = "run.config",
) -> int:
    """Builds a LazyFox Command and captures the config, stdout and -err, and time benchmarking
    Returns the status code of the LazyFox Command execution"""
    config = locals()
    config["BinaryHash"] = os.popen(f"sha1sum {binary_path}").read()

    with open(os.path.join(output_directory, config_file), "w") as f:
        json.dump(config, f)

    command = f"{binary_path} --input-graph {input_graph_path} --output-dir {output_directory} " \
              f"--queue-size {queue_size} --thread-count {thread_count}"
    if not dumping:
        command += " --disable-dumping "

    # Add logging of stderr and stdout
    command = command + f" 2>&1 > {os.path.join(output_directory, log_file)}"

    # Add time benchmark to the command
    if os.path.exists("/usr/bin/time"):
        benchmark_prefix = f"/usr/bin/time -v -o {os.path.join(output_directory, benchmark_file)}"

        command = f"{benchmark_prefix} {command}"
    else:
        print("'/usr/bin/time' not found, running without benchmark!")

    return os.system(command)


def run_on_dataset(dataset: str):
    """Converts the given dataset to a LazyFOX readable edge list,
    creates a FOX baseline, and
    creates multiple LazyFOX results with different queue sizes"""
    g = graphs[dataset]

    # Replace node labels with int labels, so lazyFox can use them
    os.makedirs(f"{OUTDIR}/edgelists", exist_ok=True)
    nx.relabel_nodes(g, mapping={node: i for i, node in enumerate(g.nodes)}, copy=False)
    nx.write_edgelist(g, f"{OUTDIR}/edgelists/{dataset}.edgelist", data=False)

    # ### LazyFox Runs ###
    run_dir = f"{OUTDIR}/{dataset}"
    os.mkdir(run_dir)

    # Baseline: FOX Run
    os.mkdir(f"{run_dir}/FOX_baseline")
    run_lazyfox(
        binary_path=LAZYFOX_BINARY,
        input_graph_path=f"{OUTDIR}/edgelists/{dataset}.edgelist",
        output_directory=f"{run_dir}/FOX_baseline",
        queue_size=1,
        thread_count=1,
        dumping=True,
    )

    # LazyFOX runs

    pbar = tqdm(QUEUE_SIZES, desc="Iterating queue sizes", leave=False)
    for queue_size in pbar:
        pbar.set_description(f"Iterating queue sizes ({queue_size})")

        os.mkdir(f"{run_dir}/LazyFOX_{queue_size}")
        run_lazyfox(
            binary_path=LAZYFOX_BINARY,
            input_graph_path=f"{OUTDIR}/edgelists/{dataset}.edgelist",
            output_directory=f"{run_dir}/LazyFOX_{queue_size}",
            queue_size=queue_size,
            thread_count=1,
            dumping=True,
        )


def main():
    pbar = tqdm(graphs, desc="Iterating graphs")
    for name in pbar:
        pbar.set_description(f"Iterating graphs ({name})")
        run_on_dataset(name)


if __name__ == "__main__":
    main()
