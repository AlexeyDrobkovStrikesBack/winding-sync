"""
Validation for angle-binned radial pitch.
Run with: python -m winding_sync.validate_pitch

The claim under test is that a single global winding-vs-radius regression is not
a reliable way to recover a scroll's radial pitch, and that the median of robust
per-wedge slopes is. All four checks run on synthetic fields with a known pitch,
so they need no bucket data and no ground-truth meshes -- a reviewer can run them
directly. (The empirical cross-scroll behaviour of the recovered field -- that it
under-resolves true density by a scroll-dependent factor -- is a property of the
CT and the solver, not of this primitive, and is documented in the module and PR,
not asserted here.)

Four checks:

1. Exactness. On a clean concentric field the returned pitch must match the
   planted pitch to within a few percent.
2. Eccentricity robustness. On an eccentric (elliptical) field a single global
   fit is confounded; median-of-wedges must stay close to the planted pitch and
   beat the global fit by a wide margin -- the reason the module exists.
3. Umbilicus refinement. The heuristic must lower per-wedge pitch dispersion and
   move toward the true centre where a concentric centre exists, and must NOT
   manufacture a low dispersion on a field with no real centre.
4. Helpers and graceful limits. The robust line fit survives gross outliers, and
   an under-populated field degrades to None rather than a fabricated number.
"""

from __future__ import annotations

import numpy as np

from .pitch import angle_binned_pitch, refine_umbilicus_pitch, theil_sen


def _circular_field(pitch_um, voxel_um, center, n_turns=8, per_turn=1500,
                    noise_w=0.03, seed=0):
    """Clean concentric spiral: w = r_um / pitch + theta / 2pi, circular."""
    rng = np.random.default_rng(seed)
    cy, cx = center
    slope = voxel_um / pitch_um                      # d(winding)/d(radius_px)
    r_max = n_turns * pitch_um / voxel_um
    n = n_turns * per_turn
    r_px = rng.uniform(0.05 * r_max, r_max, n)
    th = rng.uniform(-np.pi, np.pi, n)
    w = r_px * slope + th / (2 * np.pi) + rng.normal(0, noise_w, n)
    y = cy + r_px * np.sin(th)
    x = cx + r_px * np.cos(th)
    return np.column_stack([y, x]).astype(float), w.astype(float)


def _elliptic_field(pitch_um, voxel_um, center, ecc=2.2, n_turns=8,
                    per_turn=1500, noise_w=0.03, seed=0):
    """Eccentric concentric field: sheets are ellipses of axis ratio ``ecc``.

    Winding is linear in the elliptic radius, so along a fixed direction the
    Euclidean radius that ``angle_binned_pitch`` sees spans a wide range as a
    single winding sweeps in angle -- exactly the confound that collapses a
    global radius-vs-winding fit on real cross-sections (a Grand-Prize scroll
    measures circularity 2.2).
    """
    rng = np.random.default_rng(seed)
    cy, cx = center
    A = np.sqrt(ecc)                                 # x semi-axis scale
    B = 1.0 / np.sqrt(ecc)                           # y semi-axis scale
    slope = voxel_um / pitch_um
    rho_max = n_turns * pitch_um / voxel_um
    n = n_turns * per_turn
    rho = rng.uniform(0.05 * rho_max, rho_max, n)    # elliptic radius, px units
    th = rng.uniform(-np.pi, np.pi, n)
    w = rho * slope + th / (2 * np.pi) + rng.normal(0, noise_w, n)
    y = cy + rho * B * np.sin(th)
    x = cx + rho * A * np.cos(th)
    return np.column_stack([y, x]).astype(float), w.astype(float)


def test_exactness(verbose=True):
    vox, pitch, centre = 9.6, 360.0, (5000.0, 5000.0)
    seed_yx, w = _circular_field(pitch, vox, centre)
    res = angle_binned_pitch(seed_yx, w, centre, vox)
    err = abs(res["pitch_um"] - pitch) / pitch
    ok = res["pitch_um"] is not None and err < 0.05 and res["n_wedges"] >= 20
    if verbose:
        print("1. EXACTNESS ON A CLEAN CONCENTRIC FIELD")
        print(f"   planted pitch {pitch:.0f} um   recovered {res['pitch_um']:.1f} um"
              f"   ({err*100:.1f}% err, {res['n_wedges']} wedges, cv {res['cv']:.3f})")
        print(f"   -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_eccentricity(verbose=True):
    vox, pitch, centre = 9.6, 360.0, (5000.0, 5000.0)
    seed_yx, w = _elliptic_field(pitch, vox, centre, ecc=2.2)
    ab = angle_binned_pitch(seed_yx, w, centre, vox)
    dy = seed_yx[:, 0] - centre[0]
    dx = seed_yx[:, 1] - centre[1]
    r_px = np.hypot(dy, dx)
    gs, _ = theil_sen(r_px, w)                       # one global fit over all angles
    g_pitch = vox / abs(gs) if abs(gs) > 1e-12 else np.inf
    ab_err = abs(ab["pitch_um"] - pitch) / pitch
    g_err = abs(g_pitch - pitch) / pitch if np.isfinite(g_pitch) else np.inf
    ok = (ab["pitch_um"] is not None and ab_err < 0.15
          and g_err > 3 * ab_err and ab["n_wedges"] >= 20)
    if verbose:
        print("2. ECCENTRICITY ROBUSTNESS  (ellipse axis ratio 2.2)")
        print(f"   planted pitch {pitch:.0f} um")
        print(f"   global fit -> implied pitch {g_pitch:.0f} um   "
              f"({g_err*100:.0f}% err, confounded)")
        print(f"   angle-binned median {ab['pitch_um']:.1f} um   "
              f"({ab_err*100:.1f}% err, {ab['n_wedges']} wedges)")
        print(f"   -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_umbilicus_refine(verbose=True):
    vox, pitch, centre = 9.6, 360.0, (5000.0, 5000.0)
    seed_yx, w = _circular_field(pitch, vox, centre, seed=1)
    off = (centre[0] + 120.0, centre[1] - 120.0)
    cv0 = angle_binned_pitch(seed_yx, w, off, vox)["cv"]
    (ry, rx), cv1 = refine_umbilicus_pitch(seed_yx, w, off, vox)
    d0 = np.hypot(off[0] - centre[0], off[1] - centre[1])
    d1 = np.hypot(ry - centre[0], rx - centre[1])
    concentric_ok = cv1 < cv0 and d1 < d0

    rng = np.random.default_rng(3)                   # no real centre to find
    rseed = rng.uniform(3000, 7000, (12000, 2))
    rw = rng.normal(0, 20, 12000)
    _, cv_rand = refine_umbilicus_pitch(rseed, rw, (5000.0, 5000.0), vox)
    limit_ok = (not np.isfinite(cv_rand)) or cv_rand > 0.5

    ok = concentric_ok and limit_ok
    if verbose:
        print("3. UMBILICUS REFINEMENT  (heuristic, found/not-found signal)")
        print(f"   concentric: CV {cv0:.3f} -> {cv1:.3f}, "
              f"centre offset {d0:.0f}px -> {d1:.0f}px")
        print(f"   no-centre random field: dispersion stays high ({cv_rand:.3f})")
        print(f"   -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def test_limits(verbose=True):
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 100, 500)
    y = 2.5 * x - 7 + rng.normal(0, 3, 500)
    y[:50] = rng.uniform(-1000, 1000, 50)            # 10% gross outliers
    s, _ = theil_sen(x, y)
    ts_ok = abs(s - 2.5) < 0.15

    sparse = np.array([[5100.0, 5100.0], [5100.0, 5200.0], [5100.0, 5300.0]])
    sw = np.array([1.0, 2.0, 3.0])
    r = angle_binned_pitch(sparse, sw, (5000.0, 5000.0), 9.6)
    none_ok = r["pitch_um"] is None

    ok = ts_ok and none_ok
    if verbose:
        print("4. HELPERS AND GRACEFUL LIMITS")
        print(f"   Theil-Sen slope on a 10%-outlier line: {s:.3f} (truth 2.5)")
        print(f"   under-populated field returns pitch_um=None: {none_ok}")
        print(f"   -> {'PASS' if ok else 'FAIL'}\n")
    return ok


def main():
    print("=" * 72)
    print("angle-binned radial pitch validation")
    print("=" * 72 + "\n")
    results = [
        test_exactness(),
        test_eccentricity(),
        test_umbilicus_refine(),
        test_limits(),
    ]
    print("=" * 72)
    print(f"{sum(results)}/{len(results)} passed")
    print("=" * 72)
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
