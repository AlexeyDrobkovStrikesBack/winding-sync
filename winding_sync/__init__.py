"""winding-sync: automatic relative winding constraints from CT."""

from __future__ import annotations

# Windows consoles default to cp1252, which cannot encode many characters that
# are perfectly ordinary elsewhere (Greek letters, arrows, superscripts). A
# single such character in a print statement raises UnicodeEncodeError and kills
# the process -- which is exactly how a cosmetic table header once took down a
# whole validation run. Degrade to "?" instead of crashing.
import sys as _sys

for _stream in (_sys.stdout, _sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except Exception:
        pass

"""Scroll triage: rank Herculaneum scroll volumes by tractability, from coarse
pyramid levels, without loading full volumes."""
from .volume import VolumeSource, estimate_umbilicus
from .pitch import angle_binned_pitch, refine_umbilicus_pitch

__version__ = "0.1.0"
