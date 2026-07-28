"""Angle-binned radial pitch of a recovered winding field.

This module addresses the open problem stated in the README under "What does
not work yet": recovering an *absolute* radial pitch (micrometres of radius per
winding) from the relative winding field, without assuming a spiral.

The naive approach -- regress winding against radius over the whole slice and
invert the slope -- fails on real cross-sections, and is a plausible cause of
the reported drift (Paris 4 counts 128 -> 184 -> 61 -> 145 across resolutions).
Two things confound a global fit:

  * Eccentricity. Grand Prize sections are crushed (one measures a circularity
    ratio of 2.20). A single winding then spans nearly the full radius range as
    the angle varies, so radius and winding are only weakly correlated globally
    and the fitted slope collapses toward zero (and can flip sign).
  * The angular ramp. The winding field carries a theta/(2*pi) term; folded into
    a global radius-vs-winding fit it further destroys the correlation.

Within a *narrow angular wedge* both effects are frozen out: a wedge sees a
small, monotone radius interval per winding, so winding is locally linear in
radius. I fit each wedge with a robust (Theil-Sen) line and take the median
per-wedge slope. This reproduces the mesh radial pitch on Paris 4 where the
global fit does not.

HONEST SCOPE -- read before trusting the number as an absolute count.
`angle_binned_pitch` returns the recovered field's *own* radial pitch. Cross-
scroll validation against Grand-Prize segment meshes (4 scrolls: PHerc 0139,
0172, Paris 4, 1667) shows this recovered pitch systematically OVER-estimates
the true pitch, i.e. the recovered field UNDER-resolves the true winding
density, by a factor that varies scroll to scroll (a noisy under-count, not a
clean law):

    recovered_density / GT_density = 0.42 - 0.67
    (expanded validation: 47 GP segments across the 4 scrolls).

An earlier 4-point read suggested this factor rose monotonically with
voxels-per-turn (0.52 -> 0.82); that did NOT survive the expanded run. The 1667
"0.82" was an artifact of an off-centre umbilicus plus wrap-segment overlap:
re-measured on the merged 25-turn surface, 1667's pitch is ~236 um (on the 225 um
corpus median) with density ~0.43-0.67, sitting with the other scrolls. The
dominant systematic is umbilicus centring, not resolution.

So this is the eccentricity-robust *measurement primitive*, not a calibrated
absolute count: the recovered field under-counts winding density by ~1.5-2.4x,
so a single fixed spacing anchor (e.g. the 225 um corpus median) cannot calibrate
absolute counts. This module deliberately does not bake in a constant anchor.

Depends only on numpy.
"""

from __future__ import annotations

import numpy as np

__all__ = ["theil_sen", "angle_binned_pitch", "refine_umbilicus_pitch"]


def theil_sen(x, y, max_pairs=2000, seed=0):
    """Theil-Sen robust line fit. Returns (slope, intercept).

    Median of pairwise slopes; robust to the outliers that a winding field
    carries near sheet damage. Subsamples to ``max_pairs`` points for speed on
    large wedges (deterministic given ``seed``).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = x.size
    if n < 2:
        return 0.0, 0.0
    idx = np.arange(n)
    if n > max_pairs:
        idx = np.random.default_rng(seed).choice(n, max_pairs, replace=False)
    xi, yi = x[idx], y[idx]
    dx = xi[:, None] - xi[None, :]
    dy = yi[:, None] - yi[None, :]
    m = np.abs(dx) > 1e-9
    if not m.any():
        return 0.0, float(np.median(y))
    slope = float(np.median(dy[m] / dx[m]))
    intercept = float(np.median(y - slope * x))
    return slope, intercept


def _wedge_slopes(seed_yx, w, umbilicus_yx, voxel_um,
                  n_bins=24, min_per_bin=30, r_min_um=200.0):
    """Per-wedge |d(winding)/d(radius_px)| slopes about ``umbilicus_yx``."""
    uy, ux = umbilicus_yx
    dy = seed_yx[:, 0] - uy
    dx = seed_yx[:, 1] - ux
    r_px = np.hypot(dy, dx)
    th = np.arctan2(dy, dx)
    keep = r_px * voxel_um > r_min_um
    edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    slopes = []
    for k in range(n_bins):
        m = keep & (th >= edges[k]) & (th < edges[k + 1])
        if m.sum() < min_per_bin:
            continue
        s, _ = theil_sen(r_px[m], w[m])
        if np.isfinite(s) and abs(s) > 1e-12:
            slopes.append(abs(s))
    return np.asarray(slopes)


def angle_binned_pitch(seed_yx, w, umbilicus_yx, voxel_um,
                       n_bins=24, min_per_bin=30, r_min_um=200.0):
    """Recovered field's own radial pitch, in micrometres of radius per winding.

    Parameters
    ----------
    seed_yx : (N, 2) array of seed positions in voxels, (y, x).
    w : (N,) recovered winding values (the L1 solution; may be relative).
    umbilicus_yx : (y, x) centre in voxels.
    voxel_um : voxel size in micrometres.
    n_bins : number of angular wedges over [-pi, pi).
    min_per_bin : minimum seeds in a wedge for it to contribute.
    r_min_um : ignore seeds nearer the centre than this (the umbilicus core is
        ill-conditioned: tiny radii, large angular spread).

    Returns
    -------
    dict with:
        pitch_um        median per-wedge radial pitch (um per winding),
                        or None if fewer than 3 wedges are usable
        pitch_um_iqr    (25th, 75th) percentile of per-wedge pitch
        cv              dispersion of per-wedge pitch (std/median), a
                        centredness diagnostic -- low means concentric
        n_wedges        number of contributing wedges
        wedge_pitch_um  per-wedge pitches (for plotting / auditing)

    See the module docstring for why the returned value is the recovered field's
    OWN pitch and over-estimates the true pitch by a resolution-dependent factor.
    """
    slopes = _wedge_slopes(seed_yx, w, umbilicus_yx, voxel_um,
                           n_bins=n_bins, min_per_bin=min_per_bin,
                           r_min_um=r_min_um)
    if slopes.size < 3:
        return {"pitch_um": None, "pitch_um_iqr": None, "cv": None,
                "n_wedges": int(slopes.size), "wedge_pitch_um": slopes}
    pitch = voxel_um / slopes  # um radius per winding, per wedge
    med = float(np.median(pitch))
    return {
        "pitch_um": med,
        "pitch_um_iqr": (float(np.percentile(pitch, 25)),
                         float(np.percentile(pitch, 75))),
        "cv": float(np.std(pitch) / med) if med else None,
        "n_wedges": int(slopes.size),
        "wedge_pitch_um": pitch,
    }


def refine_umbilicus_pitch(seed_yx, w, umbilicus0_yx, voxel_um,
                           n_bins=24, min_per_bin=30, r_min_um=200.0,
                           schedule=((600, 150), (200, 50), (80, 20))):
    """Heuristic umbilicus refinement by minimising per-wedge pitch dispersion.

    A well-centred concentric section has a near-constant radial pitch in every
    direction, so the per-wedge pitch CV is lowest near the true centre. This is
    a coarse-to-fine grid search seeded at ``umbilicus0_yx`` (typically the
    field-convergence estimate). It uses the recovered field only, so it is
    complementary to the image-based ``estimate_umbilicus`` in ``volume``.

    HONEST LIMITS. This reliably *lowers* pitch dispersion where a concentric
    centre exists, but is a heuristic, not a validated centre-finder:
      * On the recovered field alone the CV minimum can sit a few mm from the
        true (mesh-geometric) centre -- the field-convergence seed was itself
        2.4-8.6 mm off on the validation scrolls, and the recovered field is
        noisiest at large radius, which biases the dispersion surface.
      * On a collapsed section with no real umbilicus (e.g. a severely crushed scroll)
        there is nothing to find and CV does not improve. Use the returned CV as
        a found/not-found signal, not the objective value alone. A dependable
        centredness check remains ground-truth concentricity (sheets that wind
        >= 3 turns), where available.

    Returns (umbilicus_yx, cv).
    """
    def cv_at(c):
        r = angle_binned_pitch(seed_yx, w, c, voxel_um, n_bins=n_bins,
                               min_per_bin=min_per_bin, r_min_um=r_min_um)
        if r["cv"] is None or r["n_wedges"] < n_bins * 0.5:
            return np.inf
        return r["cv"]

    best = (float(umbilicus0_yx[0]), float(umbilicus0_yx[1]))
    best_cv = cv_at(best)
    for rad, step in schedule:
        centre = best
        for dy in np.arange(-rad, rad + step, step):
            for dx in np.arange(-rad, rad + step, step):
                c = (centre[0] + dy, centre[1] + dx)
                v = cv_at(c)
                if v < best_cv:
                    best_cv, best = v, c
    return best, float(best_cv)
