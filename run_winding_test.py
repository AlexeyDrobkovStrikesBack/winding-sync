#!/usr/bin/env python3
"""
run_winding_test.py -- generate winding constraints on a FULL cross-section.

Why this run matters
--------------------
Constraint generation currently scores 0.827 exact agreement on an 11.2 mm crop
of PHerc0358, with a recovered winding range of exactly 50 -- matching the
physical prediction from the corpus spacing of 225 um.

But that crop is too small to answer the question that follows. Stacking extra
probes made results worse there, because nested probes from the same seed share
their counting errors: two probes measuring the same pair agreed only 7% of the
time. Genuinely independent redundancy needs different seeds pointing in
different directions, and a full cross-section is 6.5x wider than the crop in
each dimension -- about 42x the area, and far more independent geometry.

So this run tests one thing: does exact agreement rise, hold, or fall when the
redundancy is real rather than nested?

Usage:
    python run_winding_test.py
    python run_winding_test.py --scroll PHerc1447 --level 1
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors="replace")
    except Exception:
        pass

BUCKET = "s3://vesuvius-challenge-open-data"

# Grand Prize scroll volumes, plus PHerc. Paris 4 for validation. Paths read
# from a live bucket listing, not from memory.
VOLUMES = {
    "PHerc0358": ("20250821151737-9.362um-1.2m-113keV-masked.zarr", 9.362),
    "PHerc1447": ("20250521151220-8.640um-1.2m-116keV-masked.zarr", 8.640),
    "PHerc0800": ("20250521135224-8.640um-1.2m-116keV-masked.zarr", 8.640),
    "PHercParis4": ("20260411134726-2.400um-0.2m-78keV-masked.zarr", 2.400),
}

CORPUS_SPACING_UM = 225.0  # measured across all 13 scrolls, CV 0.064


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scroll", default="PHerc0358")
    ap.add_argument("--level", type=int, default=0,
                    help="0 is native resolution, which doubles the pixels per "
                         "winding (24 vs 12) and should make peak counting far "
                         "more reliable. Costs ~2.9 GB peak. Use 1 to compare "
                         "against earlier runs.")
    ap.add_argument("--slices", type=int, default=3)
    ap.add_argument("--tile", type=int, default=4096,
                    help="tile size when a slice is too large to process whole")
    ap.add_argument("--out", type=Path, default=HERE / "output" / "winding")
    args = ap.parse_args()

    if args.scroll not in VOLUMES:
        print(f"unknown scroll; known: {', '.join(sorted(VOLUMES))}")
        return 1

    from winding_sync.constraints import (build_constraints,
                                          build_constraints_tiled, consistency)
    from winding_sync.volume import VolumeSource
    from winding_sync.solver import solve_bfs_tree, solve_l1

    name, base_um = VOLUMES[args.scroll]
    url = f"{BUCKET}/{args.scroll}/volumes/{name}"
    args.out.mkdir(parents=True, exist_ok=True)

    print(f"Opening {args.scroll} at level {args.level} ...")
    src = VolumeSource.from_ome_zarr(url, level=args.level,
                                     base_voxel_um=base_um, name=args.scroll)
    vox = src.voxel_um
    print(f"  shape {src.shape}, voxel {vox:.2f} um")
    print(f"  a full slice spans {src.shape[1]*vox/1000:.1f} mm "
          f"(the crop test used 11.2 mm)\n")

    nz = src.shape[0]
    zs = np.linspace(nz * 0.35, nz * 0.65, args.slices).astype(int)
    results = []

    for z in zs:
        print(f"=== z = {z} ===")
        t0 = time.time()
        img = src.slice_at(int(z)).astype(np.float64)
        print(f"  loaded {img.shape} in {time.time()-t0:.0f}s")

        t0 = time.time()
        # Switch to tiles when a whole-slice structure tensor would not fit.
        # Six float32 arrays of the slice size is the working set; anything
        # approaching a few GB should be tiled rather than attempted.
        px = img.shape[0] * img.shape[1]
        need_gb = px * 4 * 6 / 1e9
        if need_gb > 3.0:
            print(f"  slice would need ~{need_gb:.1f} GB whole; processing in tiles")
            g, st = build_constraints_tiled(img, vox, progress=True)
        else:
            g, st = build_constraints(img, vox)
        print(f"  {st.get('n_seeds',0)} seeds, {st.get('n_edges',0)} constraints "
              f"({st.get('edges_per_node',0):.2f}/node) in {time.time()-t0:.0f}s")
        if g.n_edges < 50:
            print("  too few constraints; skipping")
            continue

        t0 = time.time()
        w_l1 = solve_l1(g)
        t_l1 = time.time() - t0
        w_bfs = solve_bfs_tree(g)

        # Robust range alongside max-min: a handful of badly assigned nodes can
        # inflate the extremes without reflecting the field as a whole.
        p1, p99 = np.percentile(w_l1, [1, 99])
        robust_range = float(p99 - p1)
        c_l1 = consistency(g, w_l1)
        c_bfs = consistency(g, w_bfs)

        # Predict from the SCROLL's extent, not the frame's. Measuring the
        # frame overstates the expected winding count by however much empty
        # space surrounds the scroll -- which made a plausible 0.84 look like
        # a failing 0.60 on the first run.
        mat = img > np.percentile(img, 60)
        my, mx = np.nonzero(mat)
        scroll_px = max(my.max() - my.min(), mx.max() - mx.min()) if mat.any() else 0
        extent_mm = img.shape[0] * vox / 1000.0
        scroll_mm = scroll_px * vox / 1000.0
        # Winding number runs 0 at the umbilicus to N at the outer edge, so
        # N = RADIUS / spacing. An earlier version used the diameter and was
        # therefore twice too large -- which made a 2.0x overcount on PHerc0358
        # read as a perfect 1.01 match. The error flattered the method for a
        # whole session.
        predicted = (scroll_mm / 2.0) * 1000.0 / CORPUS_SPACING_UM

        print(f"  L1  solved in {t_l1:.0f}s")
        print(f"    exact {c_l1['satisfied_exactly']:.3f}  "
              f"within1 {c_l1['within_one']:.3f}  range {c_l1['winding_range']}")
        print(f"  BFS exact {c_bfs['satisfied_exactly']:.3f}  "
              f"within1 {c_bfs['within_one']:.3f}  range {c_bfs['winding_range']}")
        print(f"  scroll spans {scroll_mm:.1f} mm of the {extent_mm:.1f} mm frame")
        print(f"  physical prediction: ~{predicted:.0f} windings "
              f"(radius/spacing)")
        print(f"  recovered: max-min {c_l1['winding_range']}, "
              f"p1-p99 {robust_range:.0f}")
        print(f"  graph: {st.get('n_components','?')} components, "
              f"largest holds {st.get('largest_component_fraction',0):.1%} of seeds\n")

        results.append({
            "z": int(z),
            "voxel_um": vox,
            "extent_mm": extent_mm,
            "scroll_mm": scroll_mm,
            "n_components": int(st.get("n_components", 0)),
            "largest_component_fraction": float(st.get("largest_component_fraction", 0)),
            "n_seeds": int(st.get("n_seeds", 0)),
            "n_edges": int(st.get("n_edges", 0)),
            "edges_per_node": float(st.get("edges_per_node", 0)),
            "rejected_probes": int(st.get("rejected_probes", 0)),
            "l1_exact": c_l1["satisfied_exactly"],
            "l1_within_one": c_l1["within_one"],
            "l1_range": c_l1["winding_range"],
            "l1_range_robust": robust_range,
            "bfs_exact": c_bfs["satisfied_exactly"],
            "bfs_range": c_bfs["winding_range"],
            "predicted_windings": float(predicted),
        })

    if not results:
        print("Nothing measured.")
        return 1

    ex = np.mean([r["l1_exact"] for r in results])
    epn = np.mean([r["edges_per_node"] for r in results])
    ratio = np.mean([r["l1_range"] / r["predicted_windings"] for r in results])
    ratio_r = np.mean([r["l1_range_robust"] / r["predicted_windings"] for r in results])

    lines = [
        f"{args.scroll} full-slice winding test (level {args.level})",
        "=" * 58,
        f"field size            {results[0]['extent_mm']:.1f} mm "
        f"(crop test was 11.2 mm)",
        f"edges per node        {epn:.2f}  (crop test was 1.02)",
        f"L1 exact agreement    {ex:.3f}  (crop test was 0.827)",
        f"BFS exact agreement   {np.mean([r['bfs_exact'] for r in results]):.3f}"
        f"  (crop test was 0.047)",
        f"recovered / predicted {ratio:.2f}  (max-min; 1.00 = matches 225 um)",
        f"  same, robust p1-p99  {ratio_r:.2f}",
        f"graph components      {np.mean([r['n_components'] for r in results]):.0f}"
        f"  (1 = fully connected; more means winding cannot be tied together)",
        f"largest component     "
        f"{np.mean([r['largest_component_fraction'] for r in results]):.1%} of seeds",
        "",
        "If exact agreement ROSE, independent redundancy helps and the method",
        "scales. If it FELL, the crop was flattering and the probe counting",
        "needs work before this is trustworthy at scale.",
    ]
    table = "\n".join(lines)
    print("\n" + table)

    (args.out / f"{args.scroll}_winding_test.txt").write_text(table + "\n")
    (args.out / f"{args.scroll}_winding_test.json").write_text(
        json.dumps(results, indent=2)
    )
    print(f"\nwritten to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
