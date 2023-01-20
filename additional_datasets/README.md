This addresses the reviewers critique of LazyFOX changing the results of FOX.

We want to prove with these runs, that a certain ratio of parallelization degree and node count keeps the differences ignorable.

Our initial assumption is there is a heuristic for each graph, estimating the parallelization degree (`queue_size`) at which the LazyFOX results do not show significatn differences to the FOX results.

Such a heuristic could be `queue size = 10% of node count` or far more complex ones like 
`queue_size * avg node degree >= X % of edge count`.

This folder contains datasets and code to create a FOX baseline and LazyFOX runs with various degrees of parallelization.
Together they produce the data necessary to empirically determine such a heuristic.

The datasets are all small, as the original LazyFOX experiments only showed relevant changes for the small EU-core dataset.
