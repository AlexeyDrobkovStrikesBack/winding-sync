#!/usr/bin/env python3
"""
verify_submission.py -- prove this repo works for someone who is not you.

Run this from inside the winding-sync folder before submitting:

    python verify_submission.py

Why a separate script instead of just running the tests: this package was
extracted from a larger codebase and its imports were rewired. Your own machine
already has every dependency of that larger codebase installed, so a broken
import or a missing requirement can hide on your machine and appear only for a
reviewer who starts from a clean checkout.

So this builds a FRESH virtual environment, installs ONLY what
requirements.txt declares, and runs everything inside it. That is what a
reviewer does. Nothing here touches your existing environments.

Add --online to also stream a real scroll slice from the open-data bucket.
That takes a few minutes and needs network; without it the run is offline.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
VENV = HERE / ".verify_venv"
LINE = "=" * 72

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors="replace")
    except Exception:
        pass


class Report:
    def __init__(self):
        self.rows: list[tuple[str, bool, str]] = []

    def add(self, name: str, ok: bool, detail: str = ""):
        self.rows.append((name, ok, detail))
        mark = "PASS" if ok else "FAIL"
        print(f"  [{mark}] {name}" + (f"  --  {detail}" if detail else ""), flush=True)
        return ok

    @property
    def failed(self):
        return [r for r in self.rows if not r[1]]


def venv_python() -> Path:
    return (VENV / "Scripts" / "python.exe") if sys.platform == "win32" \
        else (VENV / "bin" / "python")


def run(args, timeout=1800):
    return subprocess.run(args, capture_output=True, text=True, timeout=timeout)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--online", action="store_true",
                    help="also stream a real scroll slice from S3 (slower)")
    ap.add_argument("--keep", action="store_true",
                    help="keep the throwaway environment afterwards")
    args = ap.parse_args()

    print(LINE)
    print("winding-sync  --  clean-room verification")
    print(LINE)
    rep = Report()

    # --- 1. the repo is complete -------------------------------------
    print("\n1. Repository contents")
    required = [
        "README.md", "LICENSE", "requirements.txt", ".gitignore",
        "run_winding_test.py", "docs/img/winding_field.png",
        "winding_sync/__init__.py", "winding_sync/solver.py",
        "winding_sync/constraints.py", "winding_sync/orientation.py",
        "winding_sync/validate.py", "winding_sync/volume.py",
    ]
    missing = [f for f in required if not (HERE / f).exists()]
    rep.add("all expected files present", not missing,
            f"missing: {', '.join(missing)}" if missing else f"{len(required)} files")

    banner = HERE / "docs/img/winding_field.png"
    if banner.exists():
        rep.add("banner image is a real PNG",
                banner.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n",
                f"{banner.stat().st_size/1e6:.2f} MB")

    txt = (HERE / "README.md").read_text(encoding="utf-8", errors="replace")
    rep.add("README documents its own limitations",
            "not yet calibrated" in txt.lower() or "not calibrated" in txt.lower())
    rep.add("LICENSE is MIT", "MIT License" in
            (HERE / "LICENSE").read_text(encoding="utf-8", errors="replace"))

    # --- 2. no stale imports from the parent codebase ----------------
    print("\n2. No leftover references to the codebase this was extracted from")
    stale = []
    for py in sorted(HERE.rglob("*.py")):
        if VENV in py.parents or py.name == Path(__file__).name:
            continue
        body = py.read_text(encoding="utf-8", errors="replace")
        for bad in ("scrolltriage", "from .laminar", "from .tracing",
                    "from .winding import", "from .metrics"):
            if bad in body:
                stale.append(f"{py.name}:{bad}")
    rep.add("no stale module references", not stale,
            "; ".join(stale[:4]) if stale else "clean")

    non_ascii = []
    for py in sorted(HERE.rglob("*.py")):
        if VENV in py.parents:
            continue
        for i, line in enumerate(
                py.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if any(ord(c) > 127 for c in line):
                non_ascii.append(f"{py.name}:{i}")
    rep.add("source is ASCII-only (Windows console safe)", not non_ascii,
            "; ".join(non_ascii[:3]) if non_ascii else "clean")

    # --- 3. clean environment ----------------------------------------
    print("\n3. Building a throwaway environment (this is the real test)")
    if VENV.exists():
        shutil.rmtree(VENV)
    t0 = time.time()
    r = run([sys.executable, "-m", "venv", str(VENV)])
    if not rep.add("virtual environment created", r.returncode == 0,
                   r.stderr.strip()[:120]):
        return 1

    py = venv_python()
    run([str(py), "-m", "pip", "install", "--upgrade", "pip", "--quiet"])
    print("     installing from requirements.txt only ...", flush=True)
    r = run([str(py), "-m", "pip", "install", "-r", str(HERE / "requirements.txt")])
    if not rep.add("dependencies install from requirements.txt alone",
                   r.returncode == 0,
                   f"{time.time()-t0:.0f}s" if r.returncode == 0
                   else r.stderr.strip().splitlines()[-1][:140]):
        return 1

    r = run([str(py), "-m", "pip", "list", "--format=freeze"])
    installed = {l.split("==")[0].lower() for l in r.stdout.splitlines()}
    print("     installed:", ", ".join(sorted(installed - {"pip", "setuptools", "wheel"})[:10]))

    # --- 4. imports and tests in that environment ---------------------
    print("\n4. Import and test, inside the clean environment")
    mods = ["winding_sync", "winding_sync.solver", "winding_sync.orientation",
            "winding_sync.constraints", "winding_sync.volume",
            "winding_sync.validate"]
    r = run([str(py), "-c", "import " + ", ".join(mods) + "; print('ok')"],
            timeout=300)
    if not rep.add("every module imports", r.returncode == 0,
                   r.stderr.strip().splitlines()[-1][:140] if r.returncode else ""):
        return 1

    r = run([str(py), "-m", "winding_sync.validate"], timeout=1200)
    passed = "4/4 passed" in r.stdout
    rep.add("validation suite passes", passed,
            next((l.strip() for l in r.stdout.splitlines() if "passed" in l), "no result"))
    if not passed:
        print("\n".join("     " + l for l in r.stdout.splitlines()[-25:]))

    # --- 5. the README's headline claims are reproducible -------------
    print("\n5. Claims in the README reproduce")
    check = r'''
import numpy as np
from winding_sync.solver import make_spiral_graph, solve_l1, solve_bfs_tree, accuracy
def sweep(rate, skips=(), seeds=7):
    a=[];b=[]
    for s in range(seeds):
        g=make_spiral_graph(12,24,gross_error_rate=rate,theta_skips=skips,seed=s)
        a.append(accuracy(solve_l1(g),g.truth)); b.append(accuracy(solve_bfs_tree(g),g.truth))
    return float(np.mean(a)), float(np.mean(b))
l1_10,bfs_10 = sweep(0.10)
l1_15,_      = sweep(0.15)
_,bfs_02     = sweep(0.02)
l1_dense,_   = sweep(0.10,(2,3,4))
print(f"L1@10%={l1_10:.3f} BFS@10%={bfs_10:.3f} L1@15%={l1_15:.3f} "
      f"BFS@2%={bfs_02:.3f} L1dense@10%={l1_dense:.3f}")
# integrality: the LP must return exact integers with no rounding
g=make_spiral_graph(12,24,gross_error_rate=0.10,seed=0)
w=solve_l1(g)
print(f"integral={bool(np.all(w==np.rint(w)))}")
'''
    r = run([str(py), "-c", check], timeout=900)
    if r.returncode == 0:
        out = r.stdout.strip()
        print("     " + out.replace("\n", "\n     "))
        vals = dict(kv.split("=") for kv in out.split()[:5])
        rep.add("L1 beats BFS at 10% gross error",
                float(vals["L1@10%"]) > float(vals["BFS@10%"]) + 0.3,
                f"{vals['L1@10%']} vs {vals['BFS@10%']}")
        rep.add("L1 still >=0.90 at 15% error", float(vals["L1@15%"]) >= 0.90,
                vals["L1@15%"])
        rep.add("BFS already below 0.95 at 2% error", float(vals["BFS@2%"]) < 0.95,
                vals["BFS@2%"])
        rep.add("redundancy improves L1", float(vals["L1dense@10%"]) >= float(vals["L1@10%"]),
                f"{vals['L1@10%']} -> {vals['L1dense@10%']}")
        rep.add("LP returns exact integers (total unimodularity)",
                "integral=True" in out)
    else:
        rep.add("README claims reproduce", False,
                r.stderr.strip().splitlines()[-1][:140])

    # --- 6. optional: real data --------------------------------------
    if args.online:
        print("\n6. Real data over the network (this takes a few minutes)")
        r = run([str(py), str(HERE / "run_winding_test.py"),
                 "--scroll", "PHerc0358", "--level", "2", "--slices", "1"],
                timeout=3600)
        ok = r.returncode == 0 and "exact agreement" in r.stdout
        rep.add("streams a real scroll and solves it", ok,
                next((l.strip() for l in r.stdout.splitlines()
                      if "exact agreement" in l), r.stderr.strip()[:120]))
        if not ok:
            print("\n".join("     " + l for l in
                            (r.stdout + r.stderr).splitlines()[-15:]))
    else:
        print("\n6. Real-data check skipped (pass --online to include it)")

    # --- verdict ------------------------------------------------------
    if not args.keep and VENV.exists():
        shutil.rmtree(VENV, ignore_errors=True)

    print("\n" + LINE)
    if rep.failed:
        print(f"{len(rep.failed)} CHECK(S) FAILED  --  do not submit yet:")
        for n, _, d in rep.failed:
            print(f"   - {n}" + (f" ({d})" if d else ""))
        print(LINE)
        return 1
    print(f"ALL {len(rep.rows)} CHECKS PASSED")
    print("A reviewer starting from a clean checkout of this repo can install")
    print("and run it. Safe to submit.")
    print(LINE)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
