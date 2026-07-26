"""
Streaming access to scroll volumes, and umbilicus estimation.

Triage never loads a full volume. Grand Prize scrolls are 200-400 GB at their
published resolution and the 2.4 um slabs run to terabytes, so every read here
is a single z-slice from a coarse pyramid level. A whole scroll is triaged from
a few hundred slices totalling well under a gigabyte.

The ``VolumeSource`` abstraction accepts either an OME-Zarr group (local path or
S3 URL) or an in-memory array, so the same triage code runs against synthetic
test scrolls and against the real bucket without modification.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np



@dataclass
class VolumeSource:
    """A z-indexed source of 2-D cross-sections at a known physical scale."""

    array: object  # anything supporting array[z] -> 2-D
    voxel_um: float  # physical voxel size AT THIS PYRAMID LEVEL
    name: str = "volume"

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(self.array.shape)  # type: ignore[return-value]

    def slice_at(self, z: int) -> np.ndarray:
        return np.asarray(self.array[z], dtype=np.float32)

    @classmethod
    def from_ome_zarr(
        cls, url: str, level: int = 2, base_voxel_um: float | None = None, **kw
    ) -> "VolumeSource":
        """Open an OME-Zarr volume at a chosen pyramid level.

        Accepts three location forms:

        * ``s3://bucket/path`` -- read anonymously via s3fs. This is the path
          the Vesuvius Challenge documents for their open bucket and is the
          most reliable.
        * ``https://...``      -- read via fsspec's HTTP filesystem.
        * a local directory    -- read directly.

        ``level`` 0 is full resolution; each level is a factor-of-2 downsample,
        so the effective voxel size is ``base_voxel_um * 2**level``. Level 2 or
        3 is normally right for triage: winding spacing stays resolved at 4-8
        voxels per period while the data volume drops 64-512x.
        """
        import zarr

        zarr_major = int(str(zarr.__version__).split(".")[0])

        # The Vesuvius volumes are Zarr v2 format (the project's own viewer
        # links use the zarr2:// scheme), which both the 2.x and 3.x libraries
        # can read -- but they take stores differently. zarr 2.x wants an
        # fsspec mapper; zarr 3.x deprecated MutableMapping stores and prefers
        # a URL plus storage options. Supporting both matters because zarr 3.x
        # requires Python >= 3.11, and 3.10 is still widely installed.
        if url.startswith("s3://"):
            if zarr_major >= 3:
                root = zarr.open_group(url, mode="r", storage_options={"anon": True})
            else:
                import s3fs

                fs = s3fs.S3FileSystem(anon=True)
                root = zarr.open_group(fs.get_mapper(url), mode="r")
        elif url.startswith(("http://", "https://")):
            if zarr_major >= 3:
                root = zarr.open_group(url, mode="r")
            else:
                import fsspec

                root = zarr.open_group(fsspec.get_mapper(url), mode="r")
        else:
            root = zarr.open_group(url, mode="r")

        arr = root[str(level)]

        voxel_um = base_voxel_um
        if voxel_um is None:
            voxel_um = _voxel_um_from_metadata(root)
        if voxel_um is None:
            raise ValueError(
                "voxel size not found in OME-Zarr metadata; pass base_voxel_um"
            )
        return cls(array=arr, voxel_um=voxel_um * (2**level), name=kw.get("name", url))


def _voxel_um_from_metadata(root) -> float | None:
    """Pull the level-0 voxel size out of OME-NGFF multiscales metadata."""
    try:
        multiscales = root.attrs["multiscales"]
        ds = multiscales[0]["datasets"][0]
        for tr in ds.get("coordinateTransformations", []):
            if tr.get("type") == "scale":
                scale = tr["scale"]
                return float(scale[-1])
    except Exception:
        return None
    return None


def sample_ray(
    img: np.ndarray, centre: tuple[float, float], angle: float, r0: int, length: int
) -> np.ndarray | None:
    """Sample intensities along a radial ray, with bilinear interpolation.

    Returns None if the ray leaves the image.
    """
    cy, cx = centre
    rs = np.arange(r0, r0 + length, dtype=np.float64)
    ys = cy + rs * np.sin(angle)
    xs = cx + rs * np.cos(angle)

    h, w = img.shape
    if ys.min() < 0 or xs.min() < 0 or ys.max() >= h - 1 or xs.max() >= w - 1:
        return None

    y0 = np.floor(ys).astype(np.intp)
    x0 = np.floor(xs).astype(np.intp)
    fy = ys - y0
    fx = xs - x0

    v = (
        img[y0, x0] * (1 - fy) * (1 - fx)
        + img[y0 + 1, x0] * fy * (1 - fx)
        + img[y0, x0 + 1] * (1 - fy) * fx
        + img[y0 + 1, x0 + 1] * fy * fx
    )
    return v.astype(np.float64)


def estimate_umbilicus(
    img: np.ndarray,
    voxel_um: float,
    ref_contrast: float,
    expected_spacing_um: float = 160.0,
    search_radius_px: int = 40,
    coarse_step: int = 8,
    n_angles: int = 24,
) -> tuple[tuple[float, float], float]:
    """Locate the scroll's central axis in a cross-section.

    The umbilicus is defined here operationally: it is the centre about which
    the windings look most concentric, i.e. the point that maximises mean
    winding modulation over rays cast in all directions. An off-centre origin
    makes rays cross windings obliquely and at varying spacing, which smears
    the spectral peak and lowers modulation -- so modulation itself is the
    objective, and no separate criterion is needed.

    Returns
    -------
    (centre, best_score)
        ``centre`` is (y, x) in voxels; ``best_score`` is the mean modulation
        achieved there.

    Note: the objective is fairly flat near the optimum, because a slightly
    off-centre ray still crosses windings at almost the same spacing (the
    obliquity error goes as 1 - cos, which is second order). So the peak
    sharpness of this objective is NOT a useful measure of how well-defined a
    scroll's axis is. Axis quality is measured instead by the stability of the
    recovered centre across z, in ``triage.axis_stability``.
    """
    h, w = img.shape
    # Seed from the centroid of papyrus-bearing material.
    thresh = float(np.percentile(img, 60))
    mask = img > thresh
    if mask.sum() < 100:
        return (h / 2.0, w / 2.0), 0.0
    ys, xs = np.nonzero(mask)
    seed = (float(ys.mean()), float(xs.mean()))

    window = window_length_voxels(voxel_um, expected_spacing_um)
    r0 = max(8, int(round(1.2 * expected_spacing_um / voxel_um)))
    angles = np.linspace(0, 2 * np.pi, n_angles, endpoint=False)

    def score(centre: tuple[float, float]) -> float:
        vals = []
        for a in angles:
            prof = sample_ray(img, centre, a, r0, window)
            if prof is None:
                continue
            m = winding_modulation(prof, voxel_um, ref_contrast)
            if m.valid and m.snr >= 3.0:
                vals.append(m.wmr)
        return float(np.mean(vals)) if len(vals) >= n_angles // 3 else 0.0

    best = seed
    best_score = score(seed)

    # Two-pass coarse-to-fine. The objective is flat near the optimum, so more
    # passes buy nothing and cost a great deal: each candidate evaluation is
    # n_angles spectral estimates.
    radius = search_radius_px
    step = coarse_step
    for _ in range(2):
        offs = np.arange(-radius, radius + 1, step)
        pass_centre = best
        for dy in offs:
            for dx in offs:
                if dy == 0 and dx == 0:
                    continue
                cand = (pass_centre[0] + float(dy), pass_centre[1] + float(dx))
                s = score(cand)
                if s > best_score:
                    best_score = s
                    best = cand
        radius = max(2, int(round(radius / 4)))
        step = max(1, int(round(step / 4)))

    return best, float(best_score)
