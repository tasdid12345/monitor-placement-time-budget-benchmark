"""Monitor placement on Oregon-1 AS subnetworks: random-restart metaheuristic vs greedy-seeded incremental swap search."""
import gzip, io, json, os, random, sys, time
import numpy as np, requests, networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED, K, N_NODES, N_INST = 42, 10, 400, 5
BUDGETS, RUN_SEEDS = (0.2, 1.0, 5.0), (42, 43, 44)
URL = "https://snap.stanford.edu/data/oregon1_010331.txt.gz"
pc = time.perf_counter


def download(url):
    for i in range(5):  # retry with backoff on 429/5xx
        try:
            r = requests.get(url, timeout=60)
            if r.status_code == 429 or r.status_code >= 500:
                raise requests.HTTPError(str(r.status_code))
            r.raise_for_status()
            return r.content
        except requests.RequestException as e:
            print("download retry", i, e)
            time.sleep(2 ** i)
    sys.exit("download failed")


def make_instances(G, rng):
    """BFS-truncated connected subnetworks (induced subgraph of a BFS prefix is connected)."""
    nodes, used, out = sorted(G.nodes()), set(), []
    while len(out) < N_INST:
        s = rng.choice(nodes)
        if s in used:
            continue
        sub = [s] + [v for _, v in nx.bfs_edges(G, s)][:N_NODES - 1]
        if len(sub) < N_NODES:
            continue
        used.update(sub)
        out.append(nx.convert_node_labels_to_integers(G.subgraph(sub)))
    return out


def precompute(G):
    """One fixed BFS-tree shortest path per node pair: edge-bitset (int) and edge-id tuple."""
    n = G.number_of_nodes()
    adj = [list(G[u]) for u in range(n)]
    eid = {}
    for i, (u, v) in enumerate(G.edges()):
        eid[(u, v)] = eid[(v, u)] = i
    bits = [[0] * n for _ in range(n)]
    paths = [[()] * n for _ in range(n)]
    for s in range(n):
        bb, pp, q = {s: 0}, {s: ()}, [s]
        for p in q:
            for t in adj[p]:
                if t not in bb:
                    e = eid[(p, t)]
                    bb[t], pp[t] = bb[p] | (1 << e), pp[p] + (e,)
                    q.append(t)
                    if t > s:
                        bits[s][t] = bits[t][s] = bb[t]
                        paths[s][t] = paths[t][s] = pp[t]
    return bits, paths


class Run:
    """Wall-clock budget + anytime trajectory of the best complete solution."""
    def __init__(self, budget):
        self.t0 = pc(); self.end = self.t0 + budget; self.best = 0; self.traj = []
    def over(self):
        return pc() > self.end
    def note(self, v):
        if v > self.best:
            self.best = v; self.traj.append((pc() - self.t0, v))
    def at(self, T):
        return max([v for t, v in self.traj if t <= T], default=0)
    def t_reach(self, thr):
        return next((t for t, v in self.traj if v >= thr), None)


def baseline(n, paths, rng, run):
    """Random restarts + first-improvement swap search, full re-evaluation from path lists."""
    def full(S):
        cov = set()
        for a in range(K):
            for b in range(a + 1, K):
                cov.update(paths[S[a]][S[b]])
        return len(cov)
    while not run.over():
        S = rng.sample(range(n), K); val = full(S); run.note(val)
        improved = True
        while improved and not run.over():
            improved = False
            moves = [(i, v) for i in range(K) for v in range(n) if v not in S]
            rng.shuffle(moves)
            for i, v in moves:
                if run.over():
                    break
                T = S.copy(); T[i] = v; nv = full(T)
                if nv > val:
                    S, val, improved = T, nv, True; run.note(val)
                    break


def improved(n, bits, rng, run):
    """Greedy marginal-gain construction, best-swap search with cached bitsets, then perturb-and-repeat."""
    def value(S):
        c = 0
        for a in range(K):
            for b in range(a + 1, K):
                c |= bits[S[a]][S[b]]
        return c.bit_count()

    def local_search(S, val):
        while True:
            bv, move = val, None
            for i in range(K):
                if run.over():
                    return S, val
                rest = S[:i] + S[i + 1:]
                base = 0
                for a in range(K - 1):
                    for b in range(a + 1, K - 1):
                        base |= bits[rest[a]][rest[b]]
                rows = [bits[u] for u in rest]
                for v in range(n):  # incremental: only the 9 new pair paths are OR-ed in
                    if v in S:
                        continue
                    c = base
                    for r in rows:
                        c |= r[v]
                    cnt = c.bit_count()
                    if cnt > bv:
                        bv, move = cnt, (i, v)
            if move is None:
                return S, val
            S = S.copy(); S[move[0]] = move[1]; val = bv; run.note(val)

    # greedy: best first pair, then add the monitor with largest marginal gain
    bp = max(((bits[s][t].bit_count(), s, t) for s in range(n) for t in range(s + 1, n)))
    S, cov = [bp[1], bp[2]], bits[bp[1]][bp[2]]
    while len(S) < K:
        if run.over():
            return
        rows = [bits[u] for u in S]; bc, bvx = -1, None
        for v in range(n):
            if v in S:
                continue
            c = cov
            for r in rows:
                c |= r[v]
            if c.bit_count() > bc:
                bc, bvx = c.bit_count(), v
        for r in rows:
            cov |= r[bvx]
        S.append(bvx)
    run.note(value(S))
    cur, curv = local_search(S, value(S))
    while not run.over():  # perturbation restarts from the incumbent
        T = cur.copy()
        for i in rng.sample(range(K), 3):
            T[i] = rng.choice([v for v in range(n) if v not in T])
        T, tv = local_search(T, value(T))
        run.note(tv)
        if tv >= curv:
            cur, curv = T, tv


def main():
    rng = random.Random(SEED); np.random.seed(SEED)
    print("downloading Oregon-1 ...")
    G = nx.read_edgelist(gzip.open(io.BytesIO(download(URL))), comments="#", nodetype=int)
    G.remove_edges_from(nx.selfloop_edges(G))
    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    print("graph", G.number_of_nodes(), "nodes", G.number_of_edges(), "edges")
    insts = make_instances(G, rng)
    res = {"baseline": [], "improved": []}  # per instance: list of Run (one per seed)
    refs, ms = [], []
    for ii, g in enumerate(insts):
        t = pc(); bits, paths = precompute(g); m = g.number_of_edges(); n = g.number_of_nodes()
        print(f"instance {ii}: n={n} m={m} precompute {pc() - t:.1f}s (shared, not in budget)")
        ms.append(m)
        for name in res:
            runs = []
            for sd in RUN_SEEDS:  # one max-budget run; best-at-T read off its trajectory
                r = random.Random(sd); run = Run(max(BUDGETS))
                baseline(n, paths, r, run) if name == "baseline" else improved(n, bits, r, run)
                runs.append(run)
            res[name].append(runs)
            print("  ", name, [round(x.at(b) / m, 4) for x in runs[:1] for b in BUDGETS])
        refs.append(max(x.best for name in res for x in res[name][ii]))
    met, per_budget = {}, {}
    for name in res:
        fr = np.array([[[x.at(b) / ms[i] for x in res[name][i]] for b in BUDGETS] for i in range(len(insts))])
        per_budget[name] = fr.mean(axis=(0, 2))  # mean over instances and seeds, per budget
        met[f"{name}_score"] = float(fr.mean())
        for b, v in zip(BUDGETS, per_budget[name]):
            met[f"{name}_frac_at_{b}s"] = float(v)
        ts = [x.t_reach(0.95 * refs[i]) for i in range(len(insts)) for x in res[name][i]]
        met[f"{name}_time_to_95pct_s"] = float(np.mean([max(BUDGETS) if t is None else t for t in ts]))
        met[f"{name}_reached_95pct_frac"] = float(np.mean([t is not None for t in ts]))
    met["improvement"] = met["improved_score"] - met["baseline_score"]
    if not all(np.isfinite(v) for v in met.values()) or met["improved_score"] <= 0:
        sys.exit("no meaningful metrics")
    print(json.dumps(met, indent=2))
    os.makedirs("figures", exist_ok=True)
    x = np.arange(len(BUDGETS))
    plt.bar(x - 0.2, per_budget["baseline"], 0.4, label="baseline")
    plt.bar(x + 0.2, per_budget["improved"], 0.4, label="improved")
    plt.xticks(x, [f"{b}s" for b in BUDGETS]); plt.ylabel("mean fraction of links covered")
    plt.xlabel("wall-clock budget"); plt.legend(); plt.tight_layout(); plt.savefig("figures/coverage_by_budget.png", dpi=120)
    out = {
        "question": "Can a greedy-seeded swap local search with incremental coverage evaluation place k monitors to cover more network links than a plain randomized-restart metaheuristic under the same tight wall-clock budgets on real public network topologies?",
        "method": f"Oregon-1 AS graph cut into {N_INST} connected ~{N_NODES}-node BFS subnetworks; k={K} monitors; objective = fraction of links in the union of (fixed BFS-tree) shortest paths between monitor pairs; budgets {list(BUDGETS)} s, seeds {list(RUN_SEEDS)}.",
        "dataset": {"name": "SNAP Oregon-1 (oregon1_010331), 5 BFS subnetworks", "source_url": URL,
                    "n_rows": G.number_of_edges(), "n_train": 0, "n_test": len(insts), "same_as_paper": False},
        "metric_name": "mean fraction of links covered (over instances, seeds, budgets)",
        "higher_is_better": True,
        "baseline": "Random restarts + first-improvement swap local search with full re-evaluation of the union of pair paths.",
        "improved": "Greedy marginal-gain construction, best-swap local search with cached path-edge bitsets and incremental evaluation, then perturbation restarts until the budget ends.",
        "metrics": met,
        "notes": "No training/tuning: all 5 subnetworks are held-out test instances (n_rows = edges of the cleaned full graph; n_test = instances), so no train/test split applies and no settings were tuned. Baseline is a re-implementation (paper algorithm not available). All-pairs path precomputation is shared and excluded from budgets. Best-at-budget values are read from one 5 s anytime trajectory per seed. Time-to-95% uses the best value over all runs of both methods as reference, censored at 5 s. Subnetworks from BFS around AS hubs are small, hub-dominated and may overlap in structure. Wall-clock results depend on the runner CPU. Only 3 seeds and 5 instances.",
    }
    with open("results.json", "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
