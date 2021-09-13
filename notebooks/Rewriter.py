"""A class to rewrite the graph, ground truth and layzfox results of a specific dataset."""
from os.path import join, exists
from typing import Set
import json

from BenchmarkRun import BenchmarkRun


dataset_file_names = {
    "eu": {"graph": "email-Eu-core.txt", "gt": "email-Eu-core-department-labels.txt"},
    "dblp": {"graph": "com-dblp.ungraph.txt", "gt": "com-dblp.all.cmty.txt"},
    "lj": {"graph": "com-lj.ungraph.txt", "gt": "com-lj.all.cmty.txt"},
}

dataset_separators = {
    "eu": " ",
    "dblp": "\t",
    "lj": "\t",
}


class Rewriter(object):
    @staticmethod
    def node_labels_are_ascending(node_labels: Set[str]):
        node_labels = sorted(list(map(int, node_labels)))
        for i, node_label in enumerate(node_labels):
            if i + 1 < len(node_labels):
                if node_labels[i + 1] - 1 != node_label:
                    return False
        return True
    
    def __init__(self, dataset_directory, output_directory):
        self.dataset_directory = dataset_directory
        self.output_directory = output_directory
    
    
    def check_dataset_validity(self, dataset):
        if dataset not in dataset_file_names.keys():
            raise ValueError(f"No filenames known for dataset {dataset}")
    
    def _load_graph(self, dataset):
        """Loades the nodes and edges from a snap graph"""
        edges = set()
        node_labels = set()
        with open(join(self.dataset_directory, dataset_file_names[dataset]["graph"])) as f:
            for line in f.readlines():
                if "#" in line:
                    continue
                edge = tuple(line.strip().split(dataset_separators[dataset]))
                if len(edge) != 2:
                    raise ValueError("Bad edge!")
                edges.add(edge)
                node_labels.add(edge[0])
                node_labels.add(edge[1])
        return node_labels, edges
    
    def _load_ground_truth_communities(self, dataset):
        gt_path = dataset_file_names[dataset]["gt"]
        gt_communities = []
        # Eu has a different way of storing gt
        if dataset == "eu":
            communities = {}
            with open(join(self.dataset_directory, gt_path)) as f:
                lines = f.readlines()
                for line in lines:
                    node, com_id = line.strip().split(" ")
                    communities[com_id] = communities.get(com_id, []) + [node]

            for key in communities.keys():
                communities[key] = set(communities[key])
            gt_communities = communities.values()
        else:
            with open(join(self.dataset_directory, gt_path)) as f:
                for line in f.readlines():
                    community = line.strip().split("\t")
                    gt_communities.append(community)
        return gt_communities
        

    def rewrite_dataset(self, dataset):
        self.check_dataset_validity(dataset)
        
        # Load the graph from file
        node_labels, edges = self._load_graph(dataset)

        # nodes are not necessarily strictly ascending
        # Create mapping from node label to node index, to rewrite communities
        mapping = {}
        for i, label in enumerate(node_labels):
            mapping[label] = i

        # Write graph with remapping
        with open(join(self.output_directory, f"rewritten_{dataset}_graph.txt"), "w") as f:
            for edge in edges:
                n1, n2 = edge
                f.write(f"{mapping[n1]} {mapping[n2]}\n")
        
        # Load ground truth communities from file
        gt_communities = self._load_ground_truth_communities(dataset)
        
        # Write ground truth communities with remapping
        with open(join(self.output_directory, f"rewritten_{dataset}_gt.txt"), "w") as f:
            for community in gt_communities:
                community = list(map(lambda n_label: str(mapping[n_label]), community))
                f.write(" ".join(community))
                f.write("\n")
        
        # Write mapping
        with open(join(self.output_directory, f"node_mapping_{dataset}.json"), "w") as f:
            json.dump(mapping, f)
    
    def rewrite_lazyfox_result(self, run: BenchmarkRun, wcc_diff=0.01):
        dataset = run.dataset
        iteration = run.get_iteration_with_different_wcc_threshold(wcc_diff)
        cluster_path = join(run.run_directory, "iterations", f"{iteration}clusters.txt")
        
        if not exists(join(self.output_directory, f"node_mapping_{dataset}.json")):
            print("No node mapping for this run's dataset yet")
            print("Creating (this might take a while)")
            self.rewrite_dataset(dataset)

        # Load node mapping
        with open(join(self.output_directory, f"node_mapping_{dataset}.json"), "r") as f:
            mapping = json.load(f)
        
        # Read communities
        lazyfox_communities = []
        with open(cluster_path) as f:
            for line in f.readlines():
                community = line.strip().split("\t")
                lazyfox_communities.append(community)
        
        # Write communities with remapping
        rewritten_path = join(run.run_directory, "iterations", f"rewritten_{iteration}clusters.txt")
        with open(rewritten_path, "w") as f:
            for community in lazyfox_communities:
                community = list(map(lambda n_label: str(mapping[n_label]), community))
                f.write(" ".join(community))
                f.write("\n")