"""
Automatic relative-winding assignment by L1 integer synchronization.

The problem
-----------
Spiral fitting needs to know, for each pair of nearby surface patches, how many
windings apart they are. Getting those constraints is currently the dominant
consumer of human annotation time, and the Grand Prize tolerates only 8
documented hours of it per scroll. So this is the binding constraint on the
whole submission, not surface accuracy.

The production tool (``find_inconsistent_windings.py``) builds a BFS spanning
tree over the patch graph, propagates winding numbers along tree edges, and
uses the remaining edges to *detect* cycles whose winding delta is nonzero.
Those inconsistent loops are then rendered for a human to adjudicate one at a
time.

Two things are wrong with that as a solver:

1. **Tree propagation is maximally fragile.** A single wrong edge high in the
   tree corrupts its entire subtree, and nothing downstream can recover,
   because the tree has no redundancy by construction. The non-tree edges --
   which carry exactly the information needed to detect and repair the error --
   are used only to raise an alarm.

2. **It reports contradictions rather than resolving them.** Every detected
   loop becomes human work.

The reformulation
-----------------
Assigning integer winding numbers from noisy pairwise differences is
*synchronization over Z*, the integer-valued sibling of phase unwrapping. Given
measurements ``d_ij ~ w_i - w_j``, recover ``w``.

Minimize total weighted disagreement in L1:

    minimize  sum_ij  c_ij * | (w_i - w_j) - d_ij |

Three properties make this the right objective rather than a heuristic:

* **It is globally optimal, and exactly solvable.** Written as a linear
  program, the constraint matrix is the graph's node-arc incidence matrix
  bordered by identity blocks. Incidence matrices are totally unimodular, and
  TU is preserved under appending identity columns -- so every vertex of the LP
  polytope is integral when the measurements are. The LP relaxation therefore
  returns exact integers with no rounding and no branch-and-bound. This is the
  same structure Costantini exploited for minimum-cost-flow phase unwrapping.

* **L1 is robust to gross errors; L2 is not.** The realistic failure is not
  small noise but a patch pair being off by a whole winding. Least squares
  squares that error and smears it across the solution to reduce it; L1 pays it
  once and leaves the rest of the graph alone. This is the difference between
  degrading gracefully and failing globally.

* **It uses every edge at once.** There is no spanning tree, so no privileged
  edges whose failure is unrecoverable. Redundant cycles become error-correcting
  rather than merely error-detecting.

The scroll's spiral geometry needs no special handling. Winding number is an
unbounded integer, not a phase modulo 2*pi, so circulation around the umbilicus
emerges from the measurements themselves: a loop encircling the axis simply
accumulates +1. That makes this strictly easier than classical phase
unwrapping, which must reason modulo the wrap.

What remains genuinely hard, and is not solved here, is producing the pairwise
measurements ``d_ij`` in the first place from CT data. This module assumes they
exist and are noisy. It removes the human from *reconciling* them.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog


@dataclass
class WindingGraph:
    """A patch adjacency graph with measured relative winding offsets.

    Attributes
    ----------
    n_nodes
        Number of patches.
    edges
        (m, 2) integer array of (i, j) node pairs.
    deltas
        (m,) integer array of measured ``w_i - w_j``.
    weights
        (m,) confidence per measurement. Higher means more trusted.
    truth
        (n,) ground-truth winding numbers, when known (synthetic only).
    """

    n_nodes: int
    edges: np.ndarray
    deltas: np.ndarray
    weights: np.ndarray
    truth: np.ndarray | None = None

    @property
    def n_edges(self) -> int:
        return len(self.edges)


def solve_bfs_tree(g: WindingGraph, root: int = 0) -> np.ndarray:
    """Baseline: propagate winding numbers along a BFS spanning tree.

    This reproduces the behaviour of the current production approach. It is
    included as the thing to beat, not as a recommendation: it consults each
    node exactly once, so any erroneous edge silently contaminates everything
    beneath it and the redundant edges that could expose the error are never
    used to correct it.
    """
    from collections import deque

    adj: list[list[tuple[int, int]]] = [[] for _ in range(g.n_nodes)]
    for (i, j), d in zip(g.edges, g.deltas):
        adj[i].append((j, -int(d)))  # w_j = w_i - d
        adj[j].append((i, int(d)))  # w_i = w_j + d

    w = np.zeros(g.n_nodes, dtype=np.int64)
    seen = np.zeros(g.n_nodes, dtype=bool)
    seen[root] = True
    q = deque([root])
    while q:
        u = q.popleft()
        for v, delta in adj[u]:
            if not seen[v]:
                seen[v] = True
                w[v] = w[u] + delta
                q.append(v)
    return w


def solve_l2(g: WindingGraph) -> np.ndarray:
    """Least-squares synchronization, for comparison.

    Included to show that the robustness of the proposed method comes from the
    L1 norm specifically, not merely from using all edges at once. L2 uses every
    edge too, and still fails under gross errors, because squaring a
    whole-winding mistake makes it worth distorting the rest of the solution to
    reduce it.
    """
    m, n = g.n_edges, g.n_nodes
    rows = np.repeat(np.arange(m), 2)
    cols = g.edges.ravel()
    vals = np.tile(np.array([1.0, -1.0]), m) * np.repeat(np.sqrt(g.weights), 2)
    B = sp.csr_matrix((vals, (rows, cols)), shape=(m, n))
    rhs = g.deltas * np.sqrt(g.weights)

    # Gauge fix: pin node 0, since w is determined only up to a constant.
    B = sp.hstack([sp.csr_matrix((m, 1)), B[:, 1:]]).tocsr()
    sol = sp.linalg.lsqr(B, rhs, atol=1e-10, btol=1e-10)[0]
    sol[0] = 0.0
    return np.rint(sol).astype(np.int64)


def solve_l1(g: WindingGraph, verbose: bool = False) -> np.ndarray:
    """Globally optimal L1 winding assignment.

    Formulated as the linear program

        min  sum_e c_e (u_e + v_e)
        s.t. w_i - w_j - u_e + v_e = d_e   for each edge e = (i, j)
             u, v >= 0
             w_0 = 0                       (gauge fix)

    The constraint matrix is ``[B | -I | I]`` with ``B`` the node-arc incidence
    matrix. ``B`` is totally unimodular and TU survives appending identity
    columns, so the LP's optimal vertices are integral whenever ``d`` is --
    the returned winding numbers are exact integers, not rounded estimates.
    """
    m, n = g.n_edges, g.n_nodes

    # Incidence block over w.
    rows = np.repeat(np.arange(m), 2)
    cols = g.edges.ravel()
    vals = np.tile(np.array([1.0, -1.0]), m)
    B = sp.csr_matrix((vals, (rows, cols)), shape=(m, n))

    I = sp.identity(m, format="csr")
    A_eq = sp.hstack([B, -I, I], format="csr")
    b_eq = g.deltas.astype(float)

    # Objective: weighted L1 slack only; w itself is free.
    c = np.concatenate([np.zeros(n), g.weights, g.weights])

    bounds = [(None, None)] * n + [(0, None)] * (2 * m)
    bounds[0] = (0, 0)  # gauge fix

    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method="highs")
    if not res.success:
        raise RuntimeError(f"LP failed: {res.message}")
    if verbose:
        print(f"    LP: {n + 2*m} vars, {m} constraints, objective {res.fun:.2f}")

    w = res.x[:n]
    # Integral by total unimodularity; rint only cleans float dust.
    return np.rint(w).astype(np.int64)


def align(estimate: np.ndarray, truth: np.ndarray) -> np.ndarray:
    """Remove the global constant, which no method can determine.

    Winding number is only ever recoverable up to an additive offset, so
    comparison must quotient it out. The modal difference is used rather than
    the mean because it is unaffected by a minority of badly wrong nodes.
    """
    diff = truth - estimate
    vals, counts = np.unique(diff, return_counts=True)
    return estimate + vals[np.argmax(counts)]


def accuracy(estimate: np.ndarray, truth: np.ndarray) -> float:
    """Fraction of nodes assigned the exactly correct winding."""
    return float(np.mean(align(estimate, truth) == truth))


def make_spiral_graph(
    n_windings: int = 12,
    per_winding: int = 24,
    gross_error_rate: float = 0.0,
    missing_edge_rate: float = 0.0,
    theta_skips: tuple[int, ...] = (),
    seed: int = 0,
) -> WindingGraph:
    """Synthesise a patch graph on a spiral with known winding numbers.

    Nodes are laid out on a cylindrical grid: ``per_winding`` angular positions
    per winding, ``n_windings`` windings deep. Edges connect angular neighbours
    on the same sheet (true delta 0) and radial neighbours one winding apart
    (true delta 1), which is the connectivity a real patch graph has.

    ``theta_skips`` adds longer-range same-winding edges (connecting patches
    2, 3, ... angular steps apart), which is how graph redundancy is varied.
    Redundancy is nearly wasted on tree propagation, which consumes only n-1
    edges however many exist, but is strongly error-correcting under L1.

    Errors are injected as whole-winding offsets rather than small noise,
    because that is the realistic failure: a ray cast between two patches
    lands on the wrong sheet. Small perturbations of an integer measurement
    are not a thing that happens.
    """
    rng = np.random.default_rng(seed)
    n = n_windings * per_winding

    def idx(w, t):
        return w * per_winding + (t % per_winding)

    edges, deltas = [], []
    for w in range(n_windings):
        for t in range(per_winding):
            # Angular neighbour, same winding. Crossing the theta seam moves
            # one winding outward, which is what makes this a spiral rather
            # than a stack of rings.
            nxt = idx(w, t + 1)
            if t + 1 < per_winding:
                edges.append((idx(w, t), nxt))
                deltas.append(0)
            elif w + 1 < n_windings:
                edges.append((idx(w, t), idx(w + 1, 0)))
                deltas.append(-1)
            # Radial neighbour, one winding out.
            if w + 1 < n_windings:
                edges.append((idx(w, t), idx(w + 1, t)))
                deltas.append(-1)

    # Longer-range same-winding edges, for redundancy.
    for w in range(n_windings):
        for t in range(per_winding):
            for skip in theta_skips:
                if t + skip < per_winding:
                    edges.append((idx(w, t), idx(w, t + skip)))
                    deltas.append(0)

    edges = np.array(edges, dtype=np.int64)
    deltas = np.array(deltas, dtype=np.int64)

    truth = np.zeros(n, dtype=np.int64)
    for w in range(n_windings):
        for t in range(per_winding):
            truth[idx(w, t)] = w

    # Verify the synthetic graph is self-consistent before corrupting it.
    resid = truth[edges[:, 0]] - truth[edges[:, 1]] - deltas
    assert np.all(resid == 0), "synthetic graph is inconsistent before noise"

    if missing_edge_rate > 0:
        keep = rng.random(len(edges)) >= missing_edge_rate
        edges, deltas = edges[keep], deltas[keep]

    weights = np.ones(len(edges))

    if gross_error_rate > 0:
        n_bad = int(round(gross_error_rate * len(edges)))
        bad = rng.choice(len(edges), size=n_bad, replace=False)
        deltas = deltas.copy()
        deltas[bad] += rng.choice([-1, 1], size=n_bad)

    return WindingGraph(n, edges, deltas, weights, truth)
