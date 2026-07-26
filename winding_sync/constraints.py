"""
Automatic winding constraints from the lamina orientation field.

This closes the gap between two pieces that already work.

On one side, ``laminar.structure_tensor`` recovers the local orientation of the
papyrus laminae anywhere in a cross-section, validated on real carbonized
material: coherence 0.77-0.83 on PHerc0358 and PHerc1447. On the other,
``winding.solve_l1`` turns noisy pairwise winding offsets into globally optimal
integer winding numbers, tolerating clustered gross errors far better than the
spanning-tree propagation currently in production.

What was missing is the step between: producing those pairwise offsets from an
image without a human drawing them. That step is the whole reason winding
annotation dominates the Grand Prize's eight-hour budget.

The method
----------
The structure tensor gives, at every pixel, the direction across the laminae.
Walk a short distance along that direction and count how many laminae are
crossed: that count *is* the winding offset between the start and end points.
No centre, no spiral, no global spacing is assumed -- only that laminae are
locally distinguishable, which is exactly what the separability metric measures.

Counting is done on the image rather than by dividing distance by an assumed
spacing, because local spacing varies fourfold across a single scroll (94-412 um
on PHerc0358). Distance-based counting would inherit that error directly;
peak-counting does not.

Errors are expected and tolerated by design. A short probe that clips the edge
of a lamina, or crosses a delaminated gap, miscounts by one. That is precisely
the gross-error regime ``solve_l1`` handles, and it is why the solver was built
first: measurements this cheap can be generated in enormous redundancy, and
under L1 redundancy is error-correcting.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage, signal

from .orientation import TYPICAL_SPACING_UM, structure_tensor
from .solver import WindingGraph


@dataclass
class TracingConfig:
    """Settings for constraint generation."""

    seed_stride_um: float = 260.0     # spacing between seed points
    probe_um: float = 900.0           # how far each probe walks across laminae
    min_coherence: float = 0.35       # ignore seeds where orientation is unreliable
    max_delta: int = 6                # reject probes claiming implausible offsets
    energy_percentile: float = 40.0   # gate on gradient energy, as in laminar.py


def layer_tangent(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Unit vector ALONG the laminae.

    ``theta`` from the structure tensor is the gradient direction, which points
    across the laminae. The laminae themselves run perpendicular to it.
    """
    return -np.sin(theta), np.cos(theta)


# Physical smoothing scale applied before peak counting.
#
# NOT ANCHORED. Read the correction below before trusting any count.
#
# CORRECTION (verified against arXiv:2512.04927, full text): the figure of
# "130 windings" used throughout earlier calibration is UNSOURCED. Henderson's
# paper states a ground-truth mesh covering approximately 30 windings of a
# central subvolume of PHerc. Paris 4 (bounding box 13600x4200x3600 voxels) --
# not a whole-scroll count. No published whole-scroll winding count for that
# scroll is known to this project. Every "match against 130" claim made during
# calibration was therefore anchored to nothing, and the apparent agreement was
# meaningless. Absolute winding counts from this code remain uncalibrated.
COUNT_SMOOTH_UM = 19.0


def count_laminae(profile: np.ndarray, voxel_um: float,
                  min_spacing_um: float = 90.0) -> int:
    """How many laminae a probe crossed.

    Counted as intensity maxima, with a minimum separation derived from the
    corpus spacing so that fibre texture within one lamina is not mistaken for
    several. Returns the count of *gaps* crossed, which is what a winding offset
    measures -- n peaks means n-1 steps between them.
    """
    if profile.size < 8:
        return 0
    # Smooth by a fixed PHYSICAL distance, not a fixed number of pixels.
    #
    # A one-pixel smooth means 19.2 um of smoothing at level 3 but only 9.6 um
    # at level 2, so finer scans kept sub-lamina fibre texture that coarser ones
    # erased, and each lamina split into several counted "peaks". The same
    # Paris 4 slice returned 128 windings at 19.2 um -- matching nothing verifiable
    # and 184 at 9.6 um. The method was resolution-dependent by construction,
    # which also explains why level 0 scored worse than level 1 on PHerc0358.
    sigma_px = max(0.8, COUNT_SMOOTH_UM / voxel_um)
    smooth = ndimage.gaussian_filter1d(profile.astype(np.float64), sigma_px)
    spread = float(np.percentile(smooth, 95) - np.percentile(smooth, 5))
    if spread < 1e-6:
        return 0
    distance = max(2, int(round(min_spacing_um / voxel_um)))
    peaks, _ = signal.find_peaks(
        smooth, distance=distance, prominence=0.20 * spread
    )
    return max(0, len(peaks) - 1)


def build_constraints_tiled(
    img: np.ndarray,
    voxel_um: float,
    config: TracingConfig | None = None,
    tile: int = 4096,
    overlap_um: float = 1400.0,
    seed: int = 0,
    progress: bool = False,
) -> tuple[WindingGraph, dict]:
    """Build constraints on a large slice by processing it in tiles.

    A full-resolution Paris 4 slice is 32693 x 32693. The structure tensor needs
    roughly six arrays of that size, which is 48 GiB in float64 and still 24 GiB
    in float32 -- so a whole-slice pass simply cannot run, and trying it wastes
    twenty minutes of download before failing.

    Tiles overlap by more than one probe length so that constraints spanning a
    boundary are still formed inside some tile. Seeds are then deduplicated by
    physical position, which also stitches the per-tile graphs into one.
    """
    cfg = config or TracingConfig()
    ov = max(8, int(round(overlap_um / voxel_um)))
    step = max(64, tile - 2 * ov)
    h, w = img.shape

    all_pts: list[np.ndarray] = []
    all_edges: list[np.ndarray] = []
    all_deltas: list[np.ndarray] = []
    all_weights: list[np.ndarray] = []
    n_rejected = 0

    for y0 in range(0, h, step):
        for x0 in range(0, w, step):
            y1, x1 = min(y0 + tile, h), min(x0 + tile, w)
            if (y1 - y0) < 64 or (x1 - x0) < 64:
                continue
            sub = np.ascontiguousarray(img[y0:y1, x0:x1])
            g, st = build_constraints(sub, voxel_um, cfg, seed=seed)
            n_rejected += int(st.get("rejected_probes", 0))
            if g.n_edges == 0:
                continue
            pts = st["seed_coords"].astype(np.float64)
            pts[:, 0] += y0
            pts[:, 1] += x0
            base = sum(len(p) for p in all_pts)
            all_pts.append(pts)
            all_edges.append(g.edges + base)
            all_deltas.append(g.deltas)
            all_weights.append(g.weights)
            if progress:
                print(f"    tile ({y0},{x0}) -> {g.n_nodes} seeds, {g.n_edges} edges")
            del sub

    if not all_pts:
        return WindingGraph(0, np.zeros((0, 2), int), np.zeros(0, int),
                            np.zeros(0)), {"n_seeds": 0}

    pts = np.vstack(all_pts)
    edges = np.vstack(all_edges)
    deltas = np.concatenate(all_deltas)
    weights = np.concatenate(all_weights)

    # Merge seeds that tiles share, by PROXIMITY rather than exact coordinate.
    #
    # Ridge snapping moves each seed onto its lamina, and the same physical
    # point can land a pixel or two apart depending on which tile's crop it was
    # snapped within. Exact-coordinate matching therefore fails to merge them,
    # leaving the tiles stitched together only by accident -- which dropped the
    # largest connected component from 98% to 64%.
    from scipy.spatial import cKDTree as _KD

    merge_radius = max(1.5, 0.25 * cfg.seed_stride_um / voxel_um)
    tree_all = _KD(pts)
    groups = tree_all.query_ball_tree(tree_all, r=merge_radius)
    remap = np.full(len(pts), -1, dtype=np.int64)
    first_list: list[int] = []
    for i, grp in enumerate(groups):
        if remap[i] != -1:
            continue
        new_id = len(first_list)
        first_list.append(i)
        for j in grp:
            if remap[j] == -1:
                remap[j] = new_id
    first = np.asarray(first_list, dtype=np.int64)
    edges = remap[edges]
    keep = edges[:, 0] != edges[:, 1]
    edges, deltas, weights = edges[keep], deltas[keep], weights[keep]
    n_nodes = int(remap.max()) + 1

    import scipy.sparse as _sp
    from scipy.sparse.csgraph import connected_components as _cc

    adj = _sp.coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
                         shape=(n_nodes, n_nodes))
    n_comp, labels = _cc(adj, directed=False)
    sizes = np.bincount(labels)

    return WindingGraph(n_nodes, edges, deltas, weights), {
        "n_seeds": n_nodes,
        "n_edges": len(edges),
        "edges_per_node": len(edges) / max(1, n_nodes),
        "rejected_probes": n_rejected,
        "n_components": int(n_comp),
        "largest_component": int(sizes.max()),
        "largest_component_fraction": float(sizes.max() / n_nodes),
        "seed_coords": pts[first],
        "tiled": True,
    }


def build_constraints(
    img: np.ndarray,
    voxel_um: float,
    config: TracingConfig | None = None,
    seed: int = 0,
) -> tuple[WindingGraph, dict]:
    """Generate a winding graph from one cross-section, with no human input.

    Nodes are seed points on the laminae; edges carry the number of laminae
    crossed between them, signed by direction. The result is ready for
    ``winding.solve_l1``.

    Returns ``(graph, stats)``.
    """
    cfg = config or TracingConfig()
    theta, coh, energy = structure_tensor(img, voxel_um)
    tx, ty = layer_tangent(theta)

    gate = (energy > np.percentile(energy, cfg.energy_percentile)) & (
        coh > cfg.min_coherence
    )
    h, w = img.shape

    stride = max(3, int(round(cfg.seed_stride_um / voxel_um)))
    ys, xs = np.mgrid[stride:h - stride:stride, stride:w - stride:stride]
    ys, xs = ys.ravel(), xs.ravel()
    keep = gate[ys, xs]
    ys, xs = ys[keep], xs[keep]
    if len(ys) < 8:
        return WindingGraph(0, np.zeros((0, 2), int), np.zeros(0, int),
                            np.zeros(0)), {"n_seeds": len(ys)}

    # Snap each seed onto the nearest lamina.
    #
    # Seeds start on a regular grid, which means only about half of them land on
    # a sheet at all -- the rest sit in the gaps between sheets. That made
    # along-lamina constraints almost impossible to form: both endpoints must be
    # on the same sheet, a ~30% proposition for two independent grid points, and
    # only 20 such edges were found on an entire crop.
    #
    # Walking a short way along the gradient direction and taking the local
    # intensity maximum puts every seed on a lamina, without giving up the
    # grid's even coverage.
    snap_px = max(2, int(round(0.5 * TYPICAL_SPACING_UM / voxel_um)))
    snap_steps = np.arange(-snap_px, snap_px + 1)
    ys_s, xs_s = [], []
    for y0, x0 in zip(ys, xs):
        a = theta[y0, x0]
        yy = y0 + snap_steps * np.sin(a)
        xx = x0 + snap_steps * np.cos(a)
        ok = (yy >= 0) & (xx >= 0) & (yy < h - 1) & (xx < w - 1)
        if ok.sum() < 3:
            ys_s.append(y0); xs_s.append(x0); continue
        vals = ndimage.map_coordinates(img, [yy[ok], xx[ok]], order=1)
        best = int(np.argmax(vals))
        ys_s.append(int(round(yy[ok][best])))
        xs_s.append(int(round(xx[ok][best])))
    ys = np.asarray(ys_s, dtype=np.intp)
    xs = np.asarray(xs_s, dtype=np.intp)

    n_nodes = len(ys)
    # Index the seeds spatially so a probe can attach to the NEAREST real seed.
    #
    # The first version snapped the probe's landing point onto the seed grid and
    # looked it up by exact coordinate. But seeds only exist where the coherence
    # gate passed, so most grid positions hold no seed and the lookup simply
    # failed -- discarding roughly a third of all probes and shattering the
    # graph into 216 components on an 11 mm crop, the largest holding 30% of
    # seeds. Constraints were being satisfied within fragments while the global
    # winding field stayed undetermined.
    from scipy.spatial import cKDTree

    seed_xy = np.stack([ys, xs], axis=1).astype(float)
    tree = cKDTree(seed_xy)
    attach_radius = 1.2 * stride

    probe_px = max(4, int(round(cfg.probe_um / voxel_um)))
    steps = np.arange(0, probe_px + 1)

    edges, deltas, weights = [], [], []
    rejected = 0

    for i, (y0, x0) in enumerate(zip(ys, xs)):
        a = theta[y0, x0]
        # Probe both ways across the laminae.
        for sign in (+1, -1):
            yy = y0 + sign * steps * np.sin(a)
            xx = x0 + sign * steps * np.cos(a)
            if yy.min() < 0 or xx.min() < 0 or yy.max() >= h - 1 or xx.max() >= w - 1:
                continue
            prof = ndimage.map_coordinates(img, [yy, xx], order=1)
            n = count_laminae(prof, voxel_um)
            if n == 0 or n > cfg.max_delta:
                rejected += 1
                continue

            # Attach to the nearest real seed within a tolerance.
            ey, ex = int(round(yy[-1])), int(round(xx[-1]))
            dist, j = tree.query([float(ey), float(ex)], k=1)
            if dist > attach_radius or int(j) == i:
                rejected += 1
                continue
            j = int(j)

            # Winding increases outward along +theta by convention; the sign of
            # the probe direction carries which side we landed on.
            edges.append((i, j))
            deltas.append(-sign * n)
            # Trust the measurement more where orientation is well defined.
            weights.append(float(min(coh[y0, x0], coh[ey, ex])))

    # Edges ALONG the laminae, carrying offset zero.
    #
    # Until now only across-lamina probes were generated, giving 1.81 edges per
    # node -- far below the 4-6 that L1 needs to correct errors at the observed
    # rate, and the reason agreement stalled around 0.67 however the probes were
    # tuned. Yet following a single sheet sideways is a far EASIER measurement
    # than counting across sheets: it requires no counting at all, only that the
    # intensity stay high along the path. These constraints are both cheaper and
    # more reliable than the ones already being used, and they were simply
    # missing.
    #
    # A pair is accepted only if the vector between the seeds aligns with the
    # local lamina tangent AND the intensity along the connecting segment never
    # drops toward the inter-lamina background, which is what would happen if
    # the path had strayed onto a neighbouring sheet.
    lamina_level = float(np.percentile(img[gate], 45))
    pairs = tree.query_pairs(r=1.8 * stride, output_type="ndarray")
    for i, j in pairs:
        yi, xi = ys[i], xs[i]
        yj, xj = ys[j], xs[j]
        dy, dx = float(yj - yi), float(xj - xi)
        norm = np.hypot(dy, dx)
        if norm < 1e-6:
            continue
        tyi, txi = -np.sin(theta[yi, xi]), np.cos(theta[yi, xi])
        align = abs((dy * tyi + dx * txi) / norm)
        if align < 0.93:          # not running along the sheet
            continue
        n_steps = max(4, int(norm))
        t_ = np.linspace(0.0, 1.0, n_steps)
        seg = ndimage.map_coordinates(
            img, [yi + t_ * dy, xi + t_ * dx], order=1
        )
        if seg.min() < lamina_level:   # dipped off the sheet
            continue
        edges.append((int(i), int(j)))
        deltas.append(0)
        weights.append(float(min(coh[yi, xi], coh[yj, xj])))

    if not edges:
        return WindingGraph(n_nodes, np.zeros((0, 2), int), np.zeros(0, int),
                            np.zeros(0)), {"n_seeds": n_nodes, "rejected": rejected}

    g = WindingGraph(
        n_nodes=n_nodes,
        edges=np.asarray(edges, dtype=np.int64),
        deltas=np.asarray(deltas, dtype=np.int64),
        weights=np.asarray(weights, dtype=float),
    )
    # Connectivity. A fragmented constraint graph is fatal in a way that is
    # easy to miss: neither solver can determine relative winding between
    # components, so each gets an arbitrary offset and the global winding range
    # becomes meaningless rather than merely inaccurate. On the first full-slice
    # run BFS returned range 0 on every slice, which was this failure showing
    # itself without being named.
    import scipy.sparse as _sp
    from scipy.sparse.csgraph import connected_components as _cc

    e = np.asarray(edges, dtype=np.int64)
    adj = _sp.coo_matrix(
        (np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(n_nodes, n_nodes)
    )
    n_comp, labels = _cc(adj, directed=False)
    sizes = np.bincount(labels)

    stats = {
        "n_seeds": n_nodes,
        "n_edges": len(edges),
        "n_components": int(n_comp),
        "largest_component": int(sizes.max()),
        "largest_component_fraction": float(sizes.max() / n_nodes),
        "edges_per_node": len(edges) / max(1, n_nodes),
        "rejected_probes": rejected,
        "median_delta": float(np.median(np.abs(deltas))),
        "seed_coords": np.stack([ys, xs], axis=1),
    }
    return g, stats


def consistency(g: WindingGraph, winding: np.ndarray) -> dict:
    """How well a winding assignment explains the measurements.

    With no ground truth available on an unread scroll, this is the check that
    remains: a correct assignment should satisfy most edges exactly, because the
    measurements are integers and the truth is an integer field. A large
    residual means the constraints are contradictory -- either the probes are
    unreliable or the surface genuinely jumps.
    """
    if g.n_edges == 0:
        return {"satisfied": float("nan")}
    resid = (winding[g.edges[:, 0]] - winding[g.edges[:, 1]]) - g.deltas
    return {
        "satisfied_exactly": float(np.mean(resid == 0)),
        "within_one": float(np.mean(np.abs(resid) <= 1)),
        "mean_abs_residual": float(np.mean(np.abs(resid))),
        "winding_range": int(winding.max() - winding.min()),
    }
