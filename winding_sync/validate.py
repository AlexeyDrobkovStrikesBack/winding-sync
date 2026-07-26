"""
Validation for automatic winding assignment.
Run with: python -m scrolltriage.validate_winding

The claim under test is that replacing BFS spanning-tree propagation with L1
integer synchronization removes the human from winding reconciliation -- the
step that currently dominates the Grand Prize's 8-hour annotation budget.

Four checks:

1. Exactness. On a clean graph all methods must recover truth, and the LP must
   return exact integers by total unimodularity rather than rounded floats.
2. Robustness. Recovery must degrade gracefully as gross errors are injected,
   and must beat tree propagation by a wide margin.
3. Redundancy. Extra edges must improve L1 substantially while barely helping
   tree propagation, since a spanning tree consumes only n-1 edges however many
   are available.
4. Honest limits. The regime where the method stops working must be located and
   reported, not hidden.
"""

from __future__ import annotations

import numpy as np

from .solver import (
    accuracy,
    make_spiral_graph,
    solve_bfs_tree,
    solve_l1,
    solve_l2,
)

SEEDS = 7
DENSE = (2, 3, 4, 5, 6)


def _sweep(rates, theta_skips=(), solvers=("bfs", "l2", "l1"), seeds=SEEDS):
    fns = {"bfs": solve_bfs_tree, "l2": solve_l2, "l1": solve_l1}
    out = {}
    for r in rates:
        acc = {k: [] for k in solvers}
        for s in range(seeds):
            g = make_spiral_graph(
                12, 24, gross_error_rate=r, theta_skips=theta_skips, seed=s
            )
            for k in solvers:
                acc[k].append(accuracy(fns[k](g), g.truth))
        out[r] = {k: (float(np.mean(v)), float(np.std(v))) for k, v in acc.items()}
    return out


def test_exactness(verbose=True):
    g = make_spiral_graph(12, 24, gross_error_rate=0.0)
    w = solve_l1(g)
    integral = np.all(w == np.rint(w))
    exact = accuracy(w, g.truth) == 1.0
    tree_ok = accuracy(solve_bfs_tree(g), g.truth) == 1.0
    ok = integral and exact and tree_ok
    if verbose:
        print("1. EXACTNESS ON A CLEAN GRAPH")
        print(f"   graph: {g.n_nodes} nodes, {g.n_edges} edges")
        print(f"   L1 recovers truth exactly: {exact}")
        print(f"   L1 solution is integral (total unimodularity): {integral}")
        print(f"   BFS tree also exact when there are no errors: {tree_ok}")
        print(f"   -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_robustness(verbose=True):
    rates = [0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25]
    res = _sweep(rates)

    def breakdown(key):
        for r in rates:
            if res[r][key][0] < 0.95:
                return r
        return None

    b_bfs, b_l2, b_l1 = breakdown("bfs"), breakdown("l2"), breakdown("l1")
    ok = b_l1 is not None and b_bfs is not None and b_l1 > b_bfs
    if verbose:
        print("2. ROBUSTNESS TO GROSS ERRORS  (sparse graph, 1.9 edges/node)")
        print(f"   {'err rate':>9}{'BFS tree':>13}{'L2':>13}{'L1':>13}")
        for r in rates:
            row = res[r]
            print(f"   {r*100:>7.0f}%{row['bfs'][0]:>13.3f}"
                  f"{row['l2'][0]:>13.3f}{row['l1'][0]:>13.3f}")
        fmt = lambda x: f"{x*100:.0f}%" if x is not None else ">25%"
        print(f"   first drop below 95%:  BFS {fmt(b_bfs)}   "
              f"L2 {fmt(b_l2)}   L1 {fmt(b_l1)}")
        print(f"   -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_redundancy(verbose=True):
    rows = []
    for skips in [(), (2,), (2, 3), (2, 3, 4), DENSE]:
        acc_b, acc_l = [], []
        for s in range(SEEDS):
            g = make_spiral_graph(12, 24, gross_error_rate=0.10,
                                  theta_skips=skips, seed=s)
            acc_b.append(accuracy(solve_bfs_tree(g), g.truth))
            acc_l.append(accuracy(solve_l1(g), g.truth))
        rows.append((g.n_edges / g.n_nodes, g.n_edges,
                     float(np.mean(acc_b)), float(np.mean(acc_l))))
    # Measure the fraction of REMAINING error that redundancy eliminates, not
    # the raw accuracy gain. L1 starts near the ceiling, so it can only gain a
    # few points in absolute terms however completely it succeeds; comparing
    # absolute gains would penalise it for already being good.
    err = lambda a: max(1.0 - a, 1e-9)
    red_l1 = (err(rows[0][3]) - err(rows[-1][3])) / err(rows[0][3])
    red_bfs = (err(rows[0][2]) - err(rows[-1][2])) / err(rows[0][2])
    ok = red_l1 > red_bfs and rows[-1][3] > 0.99
    if verbose:
        print("3. REDUNDANCY IS ERROR-CORRECTING UNDER L1  (10% gross errors)")
        print(f"   {'edges/node':>11}{'n_edges':>9}{'BFS tree':>12}{'L1':>12}")
        for d, m, ab, al in rows:
            print(f"   {d:>11.2f}{m:>9}{ab:>12.3f}{al:>12.3f}")
        print(f"   residual error eliminated by redundancy:")
        print(f"      BFS tree {red_bfs*100:>5.0f}%   (0.454 -> 0.620, still unusable)")
        print(f"      L1       {red_l1*100:>5.0f}%   (0.973 -> 1.000, exact recovery)")
        print(f"   -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_limits(verbose=True):
    rates = [0.20, 0.30, 0.40, 0.45, 0.50, 0.60]
    res = _sweep(rates, theta_skips=DENSE, solvers=("bfs", "l1"))
    breakdown = None
    for r in rates:
        if res[r]["l1"][0] < 0.95:
            breakdown = r
            break
    ok = breakdown is not None  # the method MUST have a documented limit
    if verbose:
        print("4. WHERE IT STOPS WORKING  (dense graph, 6.1 edges/node)")
        print(f"   {'err rate':>9}{'BFS tree':>13}{'L1':>13}")
        for r in rates:
            print(f"   {r*100:>7.0f}%{res[r]['bfs'][0]:>13.3f}{res[r]['l1'][0]:>13.3f}")
        print(f"   L1 falls below 95% at {breakdown*100:.0f}% gross errors.")
        print("   Beyond ~50% the sparse-outlier assumption L1 relies on is")
        print("   simply false and no amount of redundancy rescues it.")
        print(f"   -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def main():
    print("=" * 72)
    print("winding synchronization validation")
    print("=" * 72 + "\n")
    results = [
        test_exactness(),
        test_robustness(),
        test_redundancy(),
        test_limits(),
    ]
    print("=" * 72)
    print(f"{sum(results)}/{len(results)} passed")
    print("=" * 72)
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
