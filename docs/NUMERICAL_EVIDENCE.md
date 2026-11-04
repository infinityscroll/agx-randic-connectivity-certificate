# Supplementary numerical evidence

These observations use floating-point arithmetic. They are not part of the exact theorem and do not certify global optimality.

The saved runs cover connected graphs at orders 4–10 and trees at orders 10–22. All 63 completed shard headers report evaluated count equal to input count, no pruning, and no mixed-order, disconnected-graph, or eigensolver-failure flags. The summed counts match the expected counts recorded in the original merge program. At order 10 the connected-graph count is 11,716,571.

| Tree order | Recorded count | Reported numerical minimizer |
|---:|---:|---|
| 10 | 106 | P10 |
| 11 | 235 | P11 |
| 12 | 551 | P12 |
| 13 | 1,301 | P13 |
| 14 | 3,159 | DC(2,2,8) |
| 15 | 7,741 | DC(2,2,9) |
| 16 | 19,320 | DC(2,2,10) |
| 17 | 48,629 | DC(2,2,11) |
| 18 | 123,867 | DC(3,3,10) |
| 19 | 317,955 | DC(3,3,11) |
| 20 | 823,065 | DC(3,3,12) |
| 21 | 2,144,505 | DC(3,3,13) |
| 22 | 5,623,756 | DC(4,4,12) |

At order 10 the saved connected-graph winner is P10 with product about 0.481037462824; the runner-up is a one-ended comet with product about 0.490743410568.

The release review recomputed all 3,077 retained candidates with a separately written graph6 decoder and NumPy eigensolver. The largest absolute disagreement in `R`, `a`, or their product was below 2.32e-14. This is agreement between implementations on retained candidates, not a rigorous error bound for every graph examined originally.

## Public replay

`numerics/selected_scan_candidates.json` contains the top two candidates for each of the 20 completed datasets, all completed shard-header fields, and SHA-256 hashes of the original shard files. The graph6 strings completely specify those 40 graphs.

```sh
python3 -m pip install numpy
python3 numerics/replay_selected_candidates.py
```

The script checks the header arithmetic and independently recomputes those 40 candidates. It does not regenerate every graph, authenticate the original run, or certify an exact global optimum. The original file labels `computed_exactly` must be understood as numerical evaluation, not exact arithmetic; the public metadata uses the original neutral `computed` field.

Generator build/version metadata is incomplete. Original shell pipelines lacked `pipefail`, so a successful final process alone would not establish enumeration completeness. Matching counts are useful checks, but do not independently prove that every isomorphism class appeared once. The partial connected-graph scan at order 11 is excluded entirely.
