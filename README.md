# Monitor placement: greedy-seeded incremental swap search vs random-restart metaheuristic

## Motivation
"A Metaheuristic Approach for Solving Monitor Placement Problem" (2023) studies how to place monitors in a communication network so that the measured paths cover as much of the network as possible. This is a hard combinatorial problem, and in practice the search often has only a short time budget. This experiment checks whether a more engineered local search helps under such tight budgets.

## Question
Can a greedy-seeded swap local search with incremental coverage evaluation place k monitors to cover more network links than a plain randomized-restart metaheuristic under the same tight wall-clock budgets (0.2 s, 1 s, 5 s) on real public network topologies?

## Dataset
The public SNAP Oregon-1 autonomous-systems graph (gzipped edge list, downloaded at run time). It is cut into 5 connected subnetworks of about 400 nodes by BFS from different random seeds. These are held-out test instances; there is no training. This is **not** the paper's own dataset.

## Baseline vs improved
- **Baseline:** a re-implemented generic metaheuristic (the paper's exact algorithm is not given). It uses random restarts and first-improvement swap local search with full re-evaluation of the covered links.
- **Improved:** greedy marginal-gain construction, then best-swap local search with cached path-edge bitsets and incremental evaluation. It then perturbs the incumbent and repeats until the budget ends.

Both choose k=10 monitors. The objective is the fraction of links in the union of shortest paths between monitor pairs. Both run under the same wall-clock budgets with 3 random seeds. Secondary metrics are the fraction covered per budget and the time to reach 95% of the best value found.

## How to run
```
pip install -r requirements.txt
python analysis.py
```

## Results
Results are written to `results.json` (and a chart to `figures/`) when the script finishes.

---
_The code in this repository was drafted with AI assistance and runs automatically on GitHub Actions; results come from the latest run._
