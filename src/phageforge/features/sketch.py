"""Tetranucleotide-frequency (TNF) genome sketches for phages.

Why TNF rather than MinHash: there are exactly ``4^4 = 256`` tetranucleotides,
which lands precisely on ``SKETCH_DIMS`` without an arbitrary projection, and
TNF is a long-established genomic signature -- compositional bias is
phylogenetically informative and stable across genomes of different length.

Both strands are counted. A FASTA record's strand orientation is arbitrary for
double-stranded phage DNA, so counting only the forward strand would make the
sketch depend on how the assembler happened to emit the contig.
"""

from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path

import numpy as np

from phageforge.config import SKETCH_DIMS

K = 4
N_KMERS = 4**K  # 256

_BASE_CODES = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    _BASE_CODES[ord(_b)] = _i
    _BASE_CODES[ord(_b.lower())] = _i

#: Complement in code space: A(0)<->T(3), C(1)<->G(2)
_COMPLEMENT = np.array([3, 2, 1, 0], dtype=np.int8)

if N_KMERS != SKETCH_DIMS:  # pragma: no cover - guards a config edit
    raise ValueError(f"TNF produces {N_KMERS} dims but SKETCH_DIMS is {SKETCH_DIMS}")


def read_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """Yield ``(header, sequence)`` pairs from a FASTA file."""
    header: str | None = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def _codes(sequence: str) -> np.ndarray:
    """Encode a sequence as base codes 0-3, with -1 for ambiguous bases."""
    raw = np.frombuffer(sequence.encode("ascii", "ignore"), dtype=np.uint8)
    return _BASE_CODES[raw]


def _count_kmers(codes: np.ndarray, counts: np.ndarray) -> None:
    """Accumulate k-mer counts for one strand into ``counts`` in place."""
    if codes.size < K:
        return
    # Sliding windows of width K, vectorised.
    windows = np.lib.stride_tricks.sliding_window_view(codes, K)
    valid = (windows >= 0).all(axis=1)
    if not valid.any():
        return
    weights = (4 ** np.arange(K - 1, -1, -1)).astype(np.int64)
    indices = windows[valid].astype(np.int64) @ weights
    counts += np.bincount(indices, minlength=N_KMERS)


def tnf_vector(sequence: str) -> np.ndarray:
    """Compute a 256-dim L2-normalised tetranucleotide frequency vector."""
    counts = np.zeros(N_KMERS, dtype=np.int64)
    codes = _codes(sequence)
    _count_kmers(codes, counts)

    # Reverse complement: reverse the array, then complement valid bases.
    rc = codes[::-1].copy()
    mask = rc >= 0
    rc[mask] = _COMPLEMENT[rc[mask]]
    _count_kmers(rc, counts)

    total = counts.sum()
    if total == 0:
        raise ValueError("sequence contains no unambiguous k-mers")

    freqs = counts.astype(np.float64) / total
    norm = np.linalg.norm(freqs)
    return (freqs / norm).astype(np.float32)


def sketch_fasta(path: Path) -> np.ndarray:
    """TNF vector for a FASTA file, concatenating all records in it."""
    sequence = "".join(seq for _, seq in read_fasta(path))
    if not sequence:
        raise ValueError(f"{path} contains no sequence")
    return tnf_vector(sequence)


def sketch_directory(directory: Path, pattern: str = "*.fna") -> dict[str, np.ndarray]:
    """Sketch every FASTA in a directory, keyed by filename stem."""
    return {path.stem: sketch_fasta(path) for path in sorted(directory.glob(pattern))}
