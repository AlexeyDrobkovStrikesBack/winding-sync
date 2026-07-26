# Release: v0.1.0

**Tag:** `v0.1.0`
**Title:** `v0.1.0 — Automatic winding constraints from CT, reconciled by L1 synchronization`

---

## Why 0.1.0 and not 1.0.0

Semantic versioning reserves `1.0.0` for a stable, defined public API. This
release has a known-uncalibrated component and an interface that will change
when that is fixed. Calling it 1.0.0 would be a claim about maturity that the
code cannot support. Under semver, `0.y.z` means exactly what is true here:
early, useful, and subject to change.

---

## Release description (paste this into the GitHub release body)

First public release.

Relative winding constraints are what stop a fitted scroll surface from
drifting onto a neighbouring wrap and rendering convincing text from the wrong
part of the scroll. Today they are drawn by hand, and contradictions between
them are resolved by hand — work that eats into the 8 documented hours of
annotation the 2027 Grand Prize allows per scroll.

This release automates both halves.

**Generation.** Constraints come straight from the CT volume. A structure
tensor gives the local orientation of the papyrus laminae everywhere, assuming
no spiral, no umbilicus, and no global winding spacing — because real scrolls
violate all three. Probes cast perpendicular to that orientation count the
laminae they cross; each count is one signed constraint. On a full PHerc0358
cross-section this yields ~20,000 seeds and ~35,000 constraints, 97% in one
connected component, with zero human input, in seconds on a laptop.

**Reconciliation.** Winding assignment from noisy pairwise offsets is
synchronization over the integers. Minimising `Σ c_ij |(w_i − w_j) − d_ij|` as
a linear program gives a constraint matrix that is totally unimodular, so the
solver returns integer solutions without rounding or branch-and-bound
(Hoffman–Kruskal 1956; the objective is Costantini's 1998, transplanted from
InSAR phase unwrapping). Verified empirically: zero integrality deviation on
the real PHerc0358 graph.

Measured against the BFS spanning-tree reconciliation in the current tooling,
accuracy first drops below 95% at **15% gross measurement error versus 2%**.
Redundancy eliminates **98%** of L1's residual error against 30% of the tree's,
because a spanning tree consumes only n−1 edges however many exist.

**A measured constant, free to anyone who needs it.** Winding spacing across
all 13 Grand Prize scrolls: **207–259 µm, median 225 µm, CV 0.064** — two
scanners, two resolutions, thirteen scrolls. Useful for probe lengths, window
sizes and sanity bounds.

### Known limitations — please read before building on this

- **Absolute winding counts are not calibrated.** The lamina counter is
  resolution-dependent; recovered counts on PHerc. Paris 4 swung 128 → 184 →
  61 → 145 across three attempted fixes. Each fix was one parameter tuned
  against one number, so tuning stopped rather than continue. Relative
  structure is sound; absolute counts are not. The proper fix is per-location
  ground truth from published segment meshes.
- **An earlier calibration target was wrong.** Development used "130 windings
  on PHerc. Paris 4" as an external anchor. That figure is unsourced —
  Henderson (WACV 2026, §3) states a ground-truth mesh covering approximately
  30 windings of a central subvolume, not a whole-scroll count. Any agreement
  previously claimed against 130 was meaningless and has been retracted
  throughout.
- **Per-constraint agreement plateaus near 0.66** on real scrolls — the same
  value across three scrolls of very different quality, including one already
  successfully read. The ceiling is in the constraint generator, not the
  specimens or the solver.
- **A trap worth knowing:** during calibration, over-smoothing cut the
  recovered count to less than half while per-constraint agreement rose to its
  best-ever 0.863. The surviving constraints were few, clean, and
  systematically wrong. Internal consistency is not correctness.

### Relationship to prior work

Henderson's diffeomorphic spiral fitting (WACV 2026) §4 contains the closest
prior approach: rays cast along local normals, a graph over surface paths,
then — in its own words — *drops nodes that are part of cycles (since these
imply a contradiction)* and majority voting. Its Limitations section names the
consequence: contradictory paths make the surface *wander between two true
windings, instead of committing to one*. This release treats a contradictory
cycle as evidence to be used rather than discarded.

### Requirements

Python 3.10–3.14, Windows/macOS/Linux. numpy, scipy, Pillow, zarr, s3fs,
fsspec. No GPU. Streams from the Vesuvius Challenge open-data bucket; nothing
loads a full volume.

### Verification

14 automated tests, including physical calibration against a closed-form
optical model and regressions pinned on real carbonized papyrus. Run
`python -m winding_sync.validate` after install. Seeds fixed
throughout.

MIT licensed.
