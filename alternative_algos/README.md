This addresses the reviewers critique of LazyFOX only being compared to FOX, not to other approaches.

This contains runs of other algorithms (`oslom v2`, `bigClam`,`core-expansion`) on our three core datasets (`EU-core`, `DBLP`, `LiveJournal`).

`oslom` and `core-expansion` struggled to compute results for the `liveJournal` dataset due to its size.
We ran `oslom` with default arg (10 runs on hierarchy 0, 50 runs on higher hierarchies), which we aborted after `68 hours`, having not finished the third run of the first ten.
We then re-ran `oslom` changing the run count to just 1 for hierarchy 0, and skipping all higher hierarchy runs.

`core-expansion` already struggled computation on `DBLP`, which took `32 hours`, which is why we did not start computation on `LiveJournal`.

# TODO
Compare the results of LazyFOX, and the other three algorihtms on the three datasets in terms of F1 and ONMI.