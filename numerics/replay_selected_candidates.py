#!/usr/bin/env python3
"""Recompute retained graph6 winners and runners-up with NumPy.

This checks 40 candidate graphs, not an exhaustive enumeration or an exact
global optimum. All eigenvalues and Randić indices use floating point.
Usage: python replay_selected_candidates.py [selected_scan_candidates.json]
"""
import json
import math
import sys
from pathlib import Path
import numpy as np


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def decode_graph6(code):
    values = [ord(c) - 63 for c in code]
    n = values[0]
    require(1 <= n <= 32, "Unsupported graph order")
    require(all(0 <= v <= 63 for v in values), "Invalid graph6 character")
    require(len(values) == 1 + (n * (n - 1) // 2 + 5) // 6, "Invalid graph6 length")
    adjacency = np.zeros((n, n), dtype=np.int64)
    bit = 0
    for j in range(1, n):
        for i in range(j):
            edge = (values[1 + bit // 6] >> (5 - bit % 6)) & 1
            adjacency[i, j] = adjacency[j, i] = edge
            bit += 1
    return adjacency


def main():
    source = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).with_name("selected_scan_candidates.json")
    data = json.loads(source.read_text())
    replay = []
    maximum_error = 0.0
    for dataset in data["datasets"]:
        total = 0
        for shard in dataset["shards"]:
            header = shard["header"]
            require(int(header["n"]) == dataset["n"], "Header order mismatch")
            require(int(header["computed"]) == int(header["count"]), "Incomplete evaluation")
            require(float(header["T"]) == 0, "Threshold pruning used")
            for field in ("mixed", "pruned_mohar", "pruned_inertia", "disconnected", "qlfail"):
                require(int(header[field]) == 0, f"Nonzero {field}")
            total += int(header["count"])
        require(total == dataset["count"], "Shard counts do not sum")
        values = []
        for rank, saved in enumerate(dataset["top_two"], 1):
            adjacency = decode_graph6(saved["graph6"])
            n = len(adjacency)
            require(n == dataset["n"], "Graph order mismatch")
            degree = adjacency.sum(axis=1)
            r = sum(1 / math.sqrt(int(degree[i] * degree[j]))
                    for i in range(n) for j in range(i) if adjacency[i, j])
            eigenvalues = np.linalg.eigvalsh(np.diag(degree) - adjacency)
            a = float(eigenvalues[1])
            require(abs(float(eigenvalues[0])) < 1e-10 and a > 1e-10, "Graph is not numerically connected")
            error = max(abs(r - saved["R"]), abs(a - saved["a"]), abs(r*a - saved["Ra"]))
            maximum_error = max(maximum_error, error)
            require(error < 1e-10, "Candidate values disagree")
            values.append({"rank":rank, "graph6":saved["graph6"], "R":r, "a":a, "Ra":r*a})
        require(values[0]["Ra"] < values[1]["Ra"], "Candidate order changed")
        replay.append({"kind":dataset["kind"], "n":dataset["n"], "candidates":values})
    print(json.dumps({"status":"PASS", "scope":"40 retained candidates and saved header arithmetic only; no global certificate",
                      "numpy_version":np.__version__, "max_absolute_disagreement":maximum_error,
                      "datasets":replay}, indent=2))


if __name__ == "__main__":
    main()
