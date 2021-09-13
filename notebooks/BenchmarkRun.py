"""Class to model a single run of LazyFox.
Used to extract runtimes and iterations numbers with varying wcc threshold"""

from os.path import join, isdir, exists
from os import listdir
import subprocess
from statistics import mean, median, stdev


dataset_stats = {
    "eu": (1_005, 25_571),
    "dblp": (317_080, 1_049_866),
    "lj": (3_997_962, 34_681_189),
}

class BenchmarkRun(object):
    def __init__(self, dataset, run_directory):
        self.base_directory = run_directory
        self.dataset = dataset
        if self.dataset not in dataset_stats.keys():
            raise ValueError(f"Missing datasets stats for '{self.dataset}'")        
        
    @property
    def run_directory(self):
        """The parent directory to where all iteration dumps are saved"""
        run_directory = [c for c in listdir(self.base_directory) if isdir(join(self.base_directory, c))][0]
        return join(self.base_directory, run_directory)
    
    @property
    def benchmark_log_lines(self):
        """The content of the file created by '/usr/bin/time -v'"""
        with open(join(self.base_directory, "bench.mark")) as f:
            return f.readlines()

    @property
    def real_timing(self):
        """Parses the benchmark file created by '/usr/bin/time -v' and returns the runtime in seconds"""
        timing_line = [l for l in self.benchmark_log_lines if "Elapsed (wall clock) time" in l][0]
        str_time = timing_line.split(" ")[-1]
        if "." in str_time:
            rest, milli = str_time.split(".")
        else:
            rest = str_time
            milli = 0
        try:
            hour, minute, second = rest.split(":")
        except ValueError:
            try:
                minute, second = rest.split(":")
                hour = 0
            except ValueError:
                raise ValueError("Unexpected amount of timing parts!")
        return int(hour) * 60 * 60 + int(minute) * 60 + int(second) + float(f"0.{milli}")
        
    @property
    def ram_peak(self):
        """IN KB, otherwise assert error"""
        ram_line = [l for l in self.benchmark_log_lines if "Maximum resident set size" in l][0]
        unit = ram_line.split("(")[1].split(")")[0]
        assert unit == "kbytes"
        return int(ram_line.split(" ")[-1].strip())
    
    @property
    def iterations(self):
        """Parses the stdout log into wcc diff and runtime per iteration"""
        with open(join(self.base_directory, "log")) as f:
            lines = f.readlines()
        
        iterations = []

        wcc_diff = None
        epoch_time = None
        for line in lines:
            if "relative change" in line:
                wcc_diff = float(line.strip().replace("relative change", ""))
            if "epoch took" in line:
                epoch_time = float(line.strip().replace("epoch took ", "").replace("s", ""))
            
            if wcc_diff is not None and epoch_time is not None:
                iteration = {
                    "wcc_diff": wcc_diff,
                    "epoch_time": epoch_time,
                }
                iterations.append(iteration)
                wcc_diff = None
                epoch_time = None

        return iterations
    
    def get_iteration_with_different_wcc_threshold(self, wcc_threshold):
        """Returns the iteration number where the given threshold was below the given one"""
        if wcc_threshold > 1 or wcc_threshold < 0.01:
            raise ValueError("Invalid WCC Threshold")
        
        iterations = self.iterations
        for i, iteration in enumerate(iterations):
            if iteration["wcc_diff"] < wcc_threshold:
                return i
    
    def get_iteration(self, num):
        """Returns wcc diff and runtime for a given iteration number"""
        return self.iteration[num]
    
    def get_iteration_timing(self, num):
        """Time passed until this iteration was done"""
        iterations = self.iterations[0:num + 1]
        return sum([iteration["epoch_time"] for iteration in iterations])
    
    def cluster_stats(self, iteration_num=None):
        """Returns some basic insights on the clusters created by LazyFox"""
        # THIS RUNS WITHOUT ANY POSTPROCESSING
        if iteration_num is None:
            iteration_num = len(self.iterations) - 1

        with open(join(self.run_directory, "iterations", f"{iteration_num}clusters.txt")) as f:
                  lines = f.readlines()
        community_count = len(lines)
        community_size_lookup = {}
        for i, line in enumerate(lines):
            community_size_lookup[i] = len(line.strip().split("\t"))
        
        lk = community_size_lookup
        
        return {
            "min": round(min(lk.values()), 2),
            "max": round(max(lk.values()), 2),
            "mean": round(mean(lk.values()), 2),
            "median": median(lk.values()),
            #"stddev": stdev(lk.values()),
            "com_count": len(lk.values()),
            "overlap": round(sum(lk.values())/ dataset_stats[self.dataset][0], 2),
        }
    
    def get_relative_epoch_time_increases(self):
        """Returns the relative increase of time spent per iteration"""
        one_off_copy = self.iterations[1:]
        diffs = []
        for i, value in enumerate(one_off_copy):
            abs_diff = value["epoch_time"] - self.iterations[i]["epoch_time"]
            diffs.append(abs_diff / self.iterations[i]["epoch_time"])
        
        return diffs
    
    def get_abs_epoch_time_increases(self):
        """Returns the absolute increase of time spent per iteration"""
        one_off_copy = self.iterations[1:]
        diffs = []
        for i, value in enumerate(one_off_copy):
            abs_diff = value["epoch_time"] - self.iterations[i]["epoch_time"]
            diffs.append(abs_diff)
        
        return diffs

    def get_saved_runtime(self, wcc_diff):
        """Returns the runtime saved by choosing a higher termination threshold"""
        last_iteration = self.get_iteration_with_different_wcc_threshold(wcc_diff)
        return sum([i["epoch_time"] for i in self.iterations[last_iteration + 1:]])
    
    def performance_stats(self):
        """Returns performance statistics of a run"""
        iterations = self.iterations
        return {
            "ram peak": self.ram_peak,
            "total time": self.real_timing,
            "first iteration runtime": self.get_iteration_timing(0),
            "iteration count": len(iterations),
            "time to epoch 4": self.get_iteration_timing(4),
            "avg time per iteration": mean([it["epoch_time"] for it in iterations]),
        }