"""Local lamina orientation from the structure tensor.

Extracted from the wider triage codebase so that winding-sync depends only on
what it actually uses. The structure tensor gives the direction across the
laminae at every point, with no assumption of a spiral, an umbilicus, or a
global winding spacing -- real scrolls violate all three.

Reference: J. Bigun and G. H. Granlund, "Optimal orientation detection of
linear symmetry", Proc. First ICCV, London, pp. 433-438, 1987.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage

TYPICAL_SPACING_UM = 225.0
GRAD_SIGMA_UM = 19.0      # gradient smoothing: below one lamina thickness
TENSOR_SIGMA_UM = 56.0    # tensor averaging: about a quarter of one spacing
ENERGY_PERCENTILE = 40.0


def structure_tensor(
    img: np.ndarray, voxel_um: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Local lamina orientation, coherence and gradient energy.

    Returns ``(theta, coherence, energy)`` where ``theta`` is the direction of
    the intensity gradient -- that is, the direction *across* the laminae, which
    is the direction to sample along.
    """
    g_sig = max(0.8, GRAD_SIGMA_UM / voxel_um)
    t_sig = max(1.5, TENSOR_SIGMA_UM / voxel_um)

    # float32 throughout. At 32693x32693 a single float64 array is 8 GiB and the
    # structure tensor needs six of them; float32 halves that, and the precision
    # is far beyond what 8-bit input data justifies.
    sm = ndimage.gaussian_filter(img.astype(np.float32), g_sig)
    gy, gx = np.gradient(sm)

    Jxx = ndimage.gaussian_filter(gx * gx, t_sig)
    Jyy = ndimage.gaussian_filter(gy * gy, t_sig)
    Jxy = ndimage.gaussian_filter(gx * gy, t_sig)

    theta = 0.5 * np.arctan2(2.0 * Jxy, Jxx - Jyy)
    trace = Jxx + Jyy
    det = Jxx * Jyy - Jxy * Jxy
    coherence = np.sqrt(np.clip(trace * trace - 4.0 * det, 0, None)) / (trace + 1e-9)
    return theta, coherence, trace
