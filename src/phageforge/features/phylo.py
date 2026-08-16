"""Bacterial strain embeddings via classical MDS on the PanACoTA distance matrix.

The Gaborieau archive ships no bacterial genomes (see ``docs/DATA_NOTES.md``), so
there is nothing to sketch with k-mers. It does ship a full strain-by-strain
phylogenetic distance matrix, which is strictly better for our purposes: rather
than approximating phylogenetic proximity with composition, we embed the real
distances directly.

Classical MDS (a.k.a. PCoA) finds coordinates whose *Euclidean* distances best
reproduce the input distances:

    B = -0.5 * J D^2 J        where J = I - (1/n) 11^T
    B = V L V^T               (eigendecomposition)
    X = V_k L_k^(1/2)         (top-k coordinates)

Because the reconstruction guarantee is Euclidean, these vectors are indexed with
``l2_norm`` similarity, not cosine.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class MDSResult:
    ids: list[str]
    coords: np.ndarray  # (n, k) float32
    eigenvalues: np.ndarray  # top-k eigenvalues
    explained: float  # fraction of positive eigenvalue mass retained
    negative_mass: float  # share of |eigenvalue| mass that was negative
    scale: float = 1.0  # uniform factor applied to coords (order-preserving)


def classical_mds(
    distances: np.ndarray,
    ids: list[str],
    n_components: int,
    *,
    scale_to_unit_rms: bool = True,
) -> MDSResult:
    """Classical MDS. ``distances`` must be square and symmetric.

    ``scale_to_unit_rms`` multiplies the coordinates by a single constant so the
    root-mean-square pairwise distance is 1. PanACoTA distances are tiny (max
    ~0.12, typical ~3e-4 between close strains), and Elasticsearch's ``l2_norm``
    score is ``1 / (1 + d^2)`` -- so unscaled, every score lands within 1e-7 of
    1.0 and float32 cannot resolve them. Scaling is a uniform similarity
    transform: it leaves the ranking identical and only restores usable numeric
    range for Stage B's similarity weighting.
    """
    d = np.asarray(distances, dtype=np.float64)
    n = d.shape[0]
    if d.shape[0] != d.shape[1]:
        raise ValueError(f"distance matrix must be square, got {d.shape}")
    if len(ids) != n:
        raise ValueError(f"{len(ids)} ids for a {n}x{n} matrix")

    # Enforce symmetry -- floating-point asymmetry in the source file would
    # otherwise produce complex eigenvalues.
    d = (d + d.T) / 2.0
    np.fill_diagonal(d, 0.0)

    # Double centering.
    squared = d**2
    j = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * j @ squared @ j
    b = (b + b.T) / 2.0

    eigenvalues, eigenvectors = np.linalg.eigh(b)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    positive = eigenvalues > 1e-9
    n_usable = int(positive.sum())
    k = min(n_components, n_usable)

    coords = eigenvectors[:, :k] * np.sqrt(eigenvalues[:k])

    # Non-Euclidean inputs yield negative eigenvalues; report how much mass that
    # is rather than discarding it silently.
    total_abs = np.abs(eigenvalues).sum()
    negative_mass = float(np.abs(eigenvalues[eigenvalues < 0]).sum() / total_abs)
    explained = float(eigenvalues[:k].sum() / eigenvalues[positive].sum())

    if k < n_components:
        # Zero-pad so every vector matches the mapping's declared dims.
        coords = np.hstack([coords, np.zeros((n, n_components - k))])

    scale = 1.0
    if scale_to_unit_rms:
        triu = np.triu_indices(n, k=1)
        rms = float(np.sqrt((d[triu] ** 2).mean()))
        if rms > 0:
            scale = 1.0 / rms
            coords = coords * scale

    return MDSResult(
        ids=list(ids),
        coords=coords.astype(np.float32),
        eigenvalues=eigenvalues[:k],
        explained=explained,
        negative_mass=negative_mass,
        scale=scale,
    )


def reconstruction_error(
    distances: np.ndarray, coords: np.ndarray, scale: float = 1.0
) -> dict[str, float]:
    """Compare embedded Euclidean distances against the originals.

    This is the quality gate: if the embedding does not reproduce the phylogeny,
    Stage B's neighbour transfer is built on sand. ``scale`` undoes the uniform
    factor applied by :func:`classical_mds` so the comparison is like for like.
    """
    d = np.asarray(distances, dtype=np.float64)
    d = (d + d.T) / 2.0

    unscaled = np.asarray(coords, dtype=np.float64) / (scale or 1.0)
    diff = unscaled[:, None, :] - unscaled[None, :, :]
    embedded = np.sqrt((diff**2).sum(axis=-1))

    triu = np.triu_indices_from(d, k=1)
    original_flat = d[triu]
    embedded_flat = embedded[triu]

    correlation = float(np.corrcoef(original_flat, embedded_flat)[0, 1])
    mae = float(np.abs(original_flat - embedded_flat).mean())
    scale = float(original_flat.mean())
    return {
        "pearson_r": correlation,
        "mae": mae,
        "mean_distance": scale,
        "relative_mae": mae / scale if scale else float("nan"),
    }
