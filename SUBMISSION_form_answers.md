# Progress Prize submission — form answers

Form: https://forms.gle/xoF5C3QsYutKP97x7
Deadline: 11:59pm Pacific, July 31st, 2026

These map onto the three Core Requirements on the prizes page. Paste into the
matching fields.

---

## Title

winding-sync — automatic relative winding constraints from CT via L1
synchronization (relative structure validated; absolute counts not yet
calibrated)

---

## Summary (short)

A tool for the winding constraints Open Problem
(scrollprize.org/open_problems/winding_annotations). It replaces both halves of
the current manual workflow: constraints are generated automatically from the
volume by probing across the local lamina orientation, and contradictions
between them are resolved globally by L1 integer synchronization instead of
being flagged for a human to adjudicate one loop at a time. The optimization is
well-behaved (totally unimodular LP, integer solutions without rounding), but
that is a statement about the solver only: the generated measurements are
directionally right and not yet calibrated in absolute count, and the output
inherits that limitation. Measured against the spanning-tree propagation in the current tooling, it
tolerates 15% gross measurement error before accuracy drops below 95%, versus
2%. On a full PHerc0358 cross-section it produces ~35,000 constraints with 97%
of seeds in one connected component, with zero human input, in seconds on a
laptop. Known limitations are documented in the README, including one
uncalibrated component and a measurement trap we hit so others don't have to.

---

## 1. Problem Identification and Solution

**The problem.** Relative winding constraints anchor scroll unrolling: they are
what stops a fitted surface drifting onto a neighbouring wrap and rendering
convincing text from the wrong part of the scroll. Today they are drawn by
hand, and the current tooling (find_inconsistent_windings.py) propagates
winding numbers along a BFS spanning tree, using remaining measurements only to
detect contradictory loops — each of which becomes annotator work. A tree is
maximally fragile (one wrong edge corrupts its whole subtree, unrecoverably)
and the redundant measurements that could repair errors are spent raising
alarms instead. This is a listed Open Problem.

**The solution — solving.** Assigning integer winding numbers from noisy
pairwise offsets is synchronization over the integers. We minimise
Σ c_ij |(w_i − w_j) − d_ij| as a linear program. The constraint matrix is the
graph's node-arc incidence matrix bordered by identity blocks — totally
unimodular, so LP vertices are integral and the solver returns integers without
rounding or branch-and-bound (the structure Costantini used for
minimum-cost-flow phase unwrapping; verified empirically on the real PHerc0358
graph, zero integrality deviation). This makes the solver a non-source of
error. It says nothing about whether the measurements fed to it are right, and
ours are uncalibrated in absolute count — see Known limitations. L1 rather than L2 because a wrong
constraint is wrong by a whole winding; least squares smears that error across
the solution, L1 pays it once.

**The solution — generating.** Constraints come from the volume with no human
input. A structure tensor gives the local lamina orientation everywhere —
assuming no spiral, no umbilicus, and no global winding spacing, since real
scrolls violate all three (circularity ratio 2.20 on one Grand Prize scroll;
local spacing varying fourfold across an 11 mm crop). Probes cast perpendicular
to that orientation count laminae crossed; each count is one signed constraint.
Individual probes are allowed to be wrong — cheap, plentiful, partly-wrong
measurements are what L1 reconciles and what tree propagation cannot survive.

**Demonstration.** Full PHerc0358 cross-section: ~20,000 seeds, ~35,000
constraints, 97% in one connected component, zero human input, seconds on a
laptop. Streaming from the open-data bucket throughout; large slices tile
automatically.

**Advantage over existing solutions, measured.** Accuracy vs injected
gross-error rate (7 seeds): first drop below 95% at BFS 2%, L2 5%, L1 15%.
Redundancy eliminates 98% of L1's residual error versus 30% of the tree's,
because a spanning tree consumes only n−1 edges however many exist. The
practical instruction inverts: collect many cheap measurements and tolerate
that a minority are wrong — which is exactly what an automatic detector
produces.

**A measured constant as a side result.** Winding spacing across all 13 Grand
Prize scrolls: 207–259 µm, median 225 µm, CV 0.064 — two scanners, two
resolutions. Useful for probe lengths, window sizes and sanity bounds; to our
knowledge unpublished.

---

## 2. Documentation

The README covers the problem, the method and the precise scope of its
mathematical properties (what the solver guarantees, and what it does not), all
measured results, install and usage for every entry point, and a prominent
"what does not work yet" section. 14 automated tests validate against ground
truth we control, including physical calibration against a closed-form optical
model and regressions pinned on real carbonized papyrus; seeds are fixed
throughout. diagnose_scroll.py produces an inspection plot (raw radial
profiles, spectrum, full-resolution crop) so anyone can examine the underlying
signal directly rather than trusting our summaries.

---

## 3. Technical Integration

- Input: OME-Zarr over s3://, https://, or local paths; zarr 2.x and 3.x;
  Python 3.10–3.14; Windows/macOS/Linux. Streaming throughout — no full volume
  is ever loaded.
- Output: a WindingGraph (nodes, edges, integer deltas, confidence weights)
  plus the L1 solution — optimal for the given measurements, which are
  themselves uncalibrated in absolute count.
- The solver is usable standalone with constraints from any source — including
  existing hand annotations — so it can replace the reconciliation step of the
  current workflow without adopting anything else here.
- Consistent formats, modular pieces, MIT licensed.

---

## Known limitations (stated so reviewers don't have to find them)

- Absolute winding counts are not calibrated. The lamina counter is
  resolution-dependent: on PHerc. Paris 4, recovered counts went
  128 → 184 → 61 → 145 across three attempted fixes. Each fix was
  one parameter tuned against one number, so we stopped rather than tune a
  fourth. The relative structure of the winding field is sound; the absolute
  count is not. The proper fix is per-location ground truth from the published
  Paris 4 segment meshes — a many-point comparison that makes overfitting
  visible. That is the next step, and help is welcome.
- Per-constraint agreement plateaus near 0.66 on real scrolls — the same value
  on three scrolls of very different quality, including one already
  successfully read — so the ceiling is in the constraint generator, not the
  specimens or the solver.
- A trap we hit, documented so others avoid it: over-smoothing cut the
  recovered count to under half of truth while per-constraint agreement rose to
  its best-ever 0.863. The surviving constraints were few, clean, and
  systematically wrong. Internal consistency is not correctness; do not
  optimise self-consistency.
- Synthetic benchmarks use regular grids; real patch graphs have clusters and
  bridges, and a bridge whose measurements are all wrong is unrecoverable
  regardless of global error rate.

---

## References

Seales, W.B., Parker, C.S., Segal, M., Tov, E., Shor, P., Porath, Y. (2016).
From damage to discovery via virtual unwrapping: Reading the scroll from
En-Gedi. Science Advances 2(9):e1601247. doi:10.1126/sciadv.1601247 — the
virtual unwrapping pipeline this field descends from.

Henderson, P. (2026). Virtually Unrolling the Herculaneum Papyri by
Diffeomorphic Spiral Fitting. WACV 2026, pp. 6401-6411. arXiv:2512.04927.
github.com/pmh47/spiral-fitting — state of the art. Sec. 4 describes the
closest prior approach to this work: rays cast along local normals, a graph of
surface paths, nodes on cycles DROPPED as contradictions, then majority voting.
Sec. 3 states the available ground truth covers ~30 windings of a central
subvolume. Sec. 7 names contradictory paths causing the surface to "wander
between two true windings" as an open limitation.

Costantini, M. (1998). A novel phase unwrapping method based on network
programming. IEEE Trans. Geoscience and Remote Sensing 36(3):813-821.
doi:10.1109/36.673674 — the L1/network-flow objective transplanted here.

Hoffman, A.J. & Kruskal, J.B. (1956). Integral Boundary Points of Convex
Polyhedra. In Linear Inequalities and Related Systems, Annals of Mathematics
Studies 38, pp. 223-246. Princeton University Press.
doi:10.1515/9781400881987-014 — total unimodularity implies integral vertices;
the theorem the solver's integrality rests on.

Singer, A. (2011). Angular synchronization by eigenvectors and semidefinite
programming. Applied and Computational Harmonic Analysis 30(1):20-36.
doi:10.1016/j.acha.2010.02.001 — the group-synchronization framing. Angular
synchronization recovers phases modulo 2pi; winding number is an unbounded
integer, so no modular ambiguity arises and the LP route is available.

Bigun, J. & Granlund, G.H. (1987). Optimal orientation detection of linear
symmetry. Proc. First ICCV, London, pp. 433-438. IEEE Computer Society — the
structure tensor, used for local lamina orientation without assuming a spiral
or umbilicus.

ScrollPrize villa repository, volume-cartographer/scripts/spiral/
find_inconsistent_windings.py — the BFS baseline measured against.
github.com/ScrollPrize/villa

Winding Constraints Open Problem:
scrollprize.org/open_problems/winding_annotations

Vesuvius Challenge open data (OME-Zarr CT volumes): scrollprize.org/data

Note on scope: these are the works this project actually stands on, verified
against primary or citing sources. Nothing is listed for appearance.

## Open source

MIT. Repository public now.
