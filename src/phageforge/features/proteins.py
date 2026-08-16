"""Recover phage receptor-binding proteins (RBPs) and embed them with ESM-2.

``RBP_list.csv`` names RBPs by Prokka locus tag (``AN17_P8_00026``) but the
archive ships **no protein sequences** (see ``docs/DATA_NOTES.md``). Recovering
them means gene-calling the phage FASTAs ourselves and taking the CDS at the
ordinal encoded in the tag.

That mapping is broadly right but not exact -- the authors' Prokka run numbers
tRNAs in the same ``locus_tag`` series and its Prodigal settings differ from
pyrodigal's meta mode, so the two gene calls drift apart at different points in
each genome, in *both* directions. A "longest protein nearby" rule would recover
~97% of tags but could just as easily grab a portal or terminase as a fibre -- a
silent wrong-protein error, which is worse than a missing one.

So the uncertainty is recorded on the document rather than averaged away, exactly
as evidence tiers are elsewhere in this project:

``high``
    Exact ordinal, and the protein is long enough to be a fibre or spike.
    Safe to embed and to use unqualified.
``heuristic``
    Exact ordinal implausibly short, but a plausible CDS sits within
    :data:`SEARCH_WINDOW` ordinals. Embedded, but flagged so the funnel can
    down-weight or exclude it.
``unresolved``
    Ordinal past the end of the genome, or nothing plausible in the window.
    Not embedded -- a missing RBP is honest, a wrong one is not.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from phageforge.config import ESM_DIMS, ESM_MODEL
from phageforge.features.sketch import read_fasta

#: Below this length a CDS cannot plausibly be a tail fibre or spike. RBP-indexed
#: proteins run to a median of ~683 aa against ~128 aa for all called CDS, so the
#: cut is well clear of the bulk of the distribution.
MIN_RBP_AA = 300

#: How far either side of the stated ordinal to look for a plausible CDS. Kept
#: deliberately tight: the wider the window, the more likely a rescue lands on an
#: unrelated structural protein.
SEARCH_WINDOW = 3

HIGH = "high"
HEURISTIC = "heuristic"
UNRESOLVED = "unresolved"


@dataclass
class CalledGene:
    """One CDS from the gene caller, numbered by its ordinal in the genome."""

    ordinal: int
    contig: str
    begin: int
    end: int
    strand: int
    protein: str

    @property
    def length(self) -> int:
        return len(self.protein.rstrip("*"))


@dataclass
class ResolvedRBP:
    """An RBP locus tag mapped onto a called CDS, with its confidence tier."""

    rbp_id: str
    phage_id: str
    rbp_type: str
    stated_ordinal: int
    confidence: str
    resolved_ordinal: int | None = None
    offset: int | None = None
    gene: CalledGene | None = None
    note: str = ""

    @property
    def protein_id(self) -> str:
        """Unique per (phage, locus tag).

        The locus tag alone is *not* unique: ``RBP_list.csv`` lists
        ``412_P4_00022`` against both ``412_P3`` and ``412_P4``. Each resolves
        against its own genome's gene call, so keying on the bare tag would let
        one phage's protein silently overwrite the other's.
        """
        return f"{self.phage_id}:{self.rbp_id}"


@dataclass
class ResolutionStats:
    total: int = 0
    by_confidence: dict[str, int] = field(default_factory=dict)
    phages_called: int = 0
    cds_called: int = 0
    rbp_lengths: list[int] = field(default_factory=list)
    all_lengths: list[int] = field(default_factory=list)

    def median_rbp_aa(self) -> float:
        return float(np.median(self.rbp_lengths)) if self.rbp_lengths else float("nan")

    def median_cds_aa(self) -> float:
        return float(np.median(self.all_lengths)) if self.all_lengths else float("nan")


# --------------------------------------------------------------- gene calling


def call_genes(path: Path) -> list[CalledGene]:
    """Gene-call one phage FASTA, numbering CDS by position across all records.

    Meta mode is used rather than per-genome training so the call matches the
    verification recorded in ``DATA_NOTES.md``; at phage genome sizes the two
    agree closely anyway.
    """
    import pyrodigal

    finder = pyrodigal.GeneFinder(meta=True)
    genes: list[CalledGene] = []
    ordinal = 0
    for header, sequence in read_fasta(path):
        for gene in finder.find_genes(sequence):
            ordinal += 1
            genes.append(
                CalledGene(
                    ordinal=ordinal,
                    contig=header,
                    begin=gene.begin,
                    end=gene.end,
                    strand=gene.strand,
                    protein=gene.translate(),
                )
            )
    return genes


def parse_locus_tag(rbp_id: str) -> tuple[str, int] | None:
    """``AN17_P8_00026`` -> ``("AN17_P8", 26)``. ``None`` if it is not a tag."""
    prefix, _, suffix = rbp_id.rpartition("_")
    if not prefix or not suffix.isdigit():
        return None
    return prefix, int(suffix)


# ------------------------------------------------------------------ resolution


def _rescue(genes: list[CalledGene], ordinal: int) -> CalledGene | None:
    """Best plausible CDS within ``SEARCH_WINDOW`` ordinals, nearest first.

    Nearest-first rather than longest-first: the failure mode we are guarding
    against is grabbing a large unrelated structural protein, and offset drift is
    small when it happens at all.
    """
    by_ordinal = {g.ordinal: g for g in genes}
    for delta in range(1, SEARCH_WINDOW + 1):
        for candidate in (ordinal - delta, ordinal + delta):
            gene = by_ordinal.get(candidate)
            if gene is not None and gene.length >= MIN_RBP_AA:
                return gene
    return None


def resolve_rbps(
    rbp_index: dict[str, list[dict[str, str]]],
    fna_dir: Path,
) -> tuple[list[ResolvedRBP], ResolutionStats]:
    """Map every RBP locus tag onto a called CDS, tiered by confidence.

    ``rbp_index`` is ``{phage_id: [{"rbp_id", "type"}]}`` as produced by
    :func:`phageforge.ingest.gaborieau._rbp_index`.
    """
    stats = ResolutionStats()
    resolved: list[ResolvedRBP] = []

    for phage_id in sorted(rbp_index):
        entries = rbp_index[phage_id]
        fasta = fna_dir / f"{phage_id}.fna"
        if not fasta.exists():
            for entry in entries:
                resolved.append(
                    ResolvedRBP(
                        rbp_id=entry["rbp_id"],
                        phage_id=phage_id,
                        rbp_type=entry.get("type", "unknown"),
                        stated_ordinal=-1,
                        confidence=UNRESOLVED,
                        note="no FASTA for this phage",
                    )
                )
            continue

        genes = call_genes(fasta)
        stats.phages_called += 1
        stats.cds_called += len(genes)
        stats.all_lengths.extend(g.length for g in genes)
        by_ordinal = {g.ordinal: g for g in genes}

        for entry in entries:
            rbp_id = entry["rbp_id"]
            parsed = parse_locus_tag(rbp_id)
            if parsed is None:
                resolved.append(
                    ResolvedRBP(
                        rbp_id=rbp_id,
                        phage_id=phage_id,
                        rbp_type=entry.get("type", "unknown"),
                        stated_ordinal=-1,
                        confidence=UNRESOLVED,
                        note="locus tag has no ordinal suffix",
                    )
                )
                continue

            tag_prefix, ordinal = parsed
            exact = by_ordinal.get(ordinal)

            # A tag whose prefix names a different phage was numbered against
            # *that* genome, so applying the ordinal here is an extra assumption
            # on top of the gene-caller drift. Never `high`, whatever it hits.
            foreign = tag_prefix != phage_id
            foreign_note = f"locus tag belongs to {tag_prefix}, resolved against {phage_id}"

            if exact is not None and exact.length >= MIN_RBP_AA:
                record = ResolvedRBP(
                    rbp_id=rbp_id,
                    phage_id=phage_id,
                    rbp_type=entry.get("type", "unknown"),
                    stated_ordinal=ordinal,
                    confidence=HEURISTIC if foreign else HIGH,
                    resolved_ordinal=ordinal,
                    offset=0,
                    gene=exact,
                    note=foreign_note if foreign else "",
                )
            else:
                rescued = _rescue(genes, ordinal)
                if rescued is not None:
                    record = ResolvedRBP(
                        rbp_id=rbp_id,
                        phage_id=phage_id,
                        rbp_type=entry.get("type", "unknown"),
                        stated_ordinal=ordinal,
                        confidence=HEURISTIC,
                        resolved_ordinal=rescued.ordinal,
                        offset=rescued.ordinal - ordinal,
                        gene=rescued,
                        note=(
                            f"exact ordinal {'absent' if exact is None else f'{exact.length} aa'}; "
                            f"recovered at offset {rescued.ordinal - ordinal:+d}"
                        ),
                    )
                else:
                    record = ResolvedRBP(
                        rbp_id=rbp_id,
                        phage_id=phage_id,
                        rbp_type=entry.get("type", "unknown"),
                        stated_ordinal=ordinal,
                        confidence=UNRESOLVED,
                        note=(
                            "ordinal past the end of the genome"
                            if exact is None
                            else f"exact ordinal is {exact.length} aa and nothing plausible within "
                            f"+/-{SEARCH_WINDOW}"
                        ),
                    )

            if record.gene is not None:
                stats.rbp_lengths.append(record.gene.length)
            resolved.append(record)

    stats.total = len(resolved)
    for record in resolved:
        stats.by_confidence[record.confidence] = stats.by_confidence.get(record.confidence, 0) + 1
    return resolved, stats


# ------------------------------------------------------------------ embedding


def _batches(items: list, size: int) -> Iterator[list]:
    for start in range(0, len(items), size):
        yield items[start : start + size]


def _cache_path() -> Path:
    from phageforge.config import CACHE_DIR

    return CACHE_DIR / f"esm-{ESM_MODEL}.npz"


def _cache_key(sequence: str) -> str:
    """Cache on the sequence itself, not the protein ID.

    The gene caller can move which CDS a tag resolves to between runs; keying on
    the sequence means a re-resolution reuses an embedding only when it really is
    the same protein.
    """
    import hashlib

    return hashlib.sha1(sequence.encode()).hexdigest()


def embed_proteins(
    sequences: dict[str, str],
    *,
    batch_size: int = 4,
    progress: bool = True,
    use_cache: bool = True,
) -> dict[str, np.ndarray]:
    """Mean-pooled ESM-2 embeddings, keyed like ``sequences``.

    Mean pooling over residue representations (excluding the BOS/EOS tokens) is
    the standard sequence-level summary for ESM-2 and is what the RBP-similarity
    literature uses. Vectors are L2-normalised because the index scores them with
    cosine.

    Runs on CPU: ``esm2_t12_35M_UR50D`` takes ~3 minutes for ~130 fibre-length
    proteins on this host, so results are cached by sequence hash under
    ``data/derived/cache``. Embedding is strictly an offline step -- never call
    it on the request path.
    """
    cache: dict[str, np.ndarray] = {}
    path = _cache_path()
    if use_cache and path.exists():
        with np.load(path) as stored:
            cache = {k: stored[k] for k in stored.files}

    out: dict[str, np.ndarray] = {}
    pending: dict[str, str] = {}
    for name, seq in sequences.items():
        hit = cache.get(_cache_key(seq))
        if hit is not None:
            out[name] = hit
        else:
            pending[name] = seq

    if progress and out:
        print(f"    {len(out)} of {len(sequences)} embeddings served from cache")
    if not pending:
        return out

    out.update(_embed_uncached(pending, batch_size=batch_size, progress=progress))

    if use_cache:
        cache.update({_cache_key(sequences[name]): out[name] for name in pending})
        np.savez_compressed(path, **cache)
    return out


def _embed_uncached(
    sequences: dict[str, str],
    *,
    batch_size: int,
    progress: bool,
) -> dict[str, np.ndarray]:
    import esm
    import torch

    loader = getattr(esm.pretrained, ESM_MODEL)
    model, alphabet = loader()
    model.eval()
    layer = model.num_layers
    if model.embed_dim != ESM_DIMS:
        raise ValueError(
            f"{ESM_MODEL} emits {model.embed_dim} dims but ESM_DIMS is {ESM_DIMS}; "
            "the pf-proteins mapping would reject these vectors"
        )

    converter = alphabet.get_batch_converter()
    # Longest first: batches are padded to their longest member, so grouping
    # similar lengths together keeps the wasted compute down.
    items = sorted(sequences.items(), key=lambda kv: -len(kv[1]))
    out: dict[str, np.ndarray] = {}
    if progress:
        print(f"    embedding {len(items)} proteins")

    with torch.no_grad():
        for n, batch in enumerate(_batches(items, batch_size), start=1):
            _, _, tokens = converter([(name, seq) for name, seq in batch])
            result = model(tokens, repr_layers=[layer], return_contacts=False)
            representations = result["representations"][layer]
            for i, (name, seq) in enumerate(batch):
                # Token 0 is BOS and token len+1 is EOS; pooling over them would
                # blend a constant into every vector.
                pooled = representations[i, 1 : len(seq) + 1].mean(0)
                vector = pooled.numpy().astype(np.float32)
                norm = float(np.linalg.norm(vector))
                out[name] = vector / norm if norm else vector
            if progress:
                print(f"    embedded {len(out)}/{len(items)}", end="\r", flush=True)
    if progress:
        print(" " * 40, end="\r")
    return out


#: Where the corpus mean direction is persisted. Any vector embedded later --
#: a newly uploaded phage at query time -- must be centred against *this* mean,
#: not against its own batch, or it lands in a different space from the index.
MEAN_VECTOR_FILE = "esm-mean.npy"


def center_vectors(vectors: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Remove the corpus mean direction, then re-normalise.

    ESM-2 embeddings are strongly **anisotropic**: every vector carries a large
    shared component, so raw cosine similarity between any two proteins in this
    corpus lands between 0.88 and 1.00 with a median of 0.97. Elasticsearch will
    happily rank on that, but it is ranking on a constant plus a whisper.

    Subtracting the mean and re-normalising restores the range to -0.75..1.00
    (median -0.06) and makes the ordering biologically legible -- two
    Felixounavirus fibres stay at 0.98 while a Tequatrovirus fibre drops from
    0.98 to 0.39. Both the Stage A RBP arm and the cocktail's receptor-overlap
    penalty are meaningless without it.

    This is the standard "all-but-the-top" correction for contextual embeddings,
    and it is a similarity transform on the whole corpus: it changes what the
    numbers can distinguish, not which protein is which.
    """
    names = list(vectors)
    matrix = np.stack([vectors[n] for n in names]).astype(np.float64)
    mean = matrix.mean(axis=0)
    centred = matrix - mean
    norms = np.linalg.norm(centred, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    centred = centred / norms
    return (
        {n: centred[i].astype(np.float32) for i, n in enumerate(names)},
        mean.astype(np.float32),
    )


def save_mean_vector(mean: np.ndarray) -> Path:
    from phageforge.config import CACHE_DIR

    path = CACHE_DIR / MEAN_VECTOR_FILE
    np.save(path, mean)
    return path


def load_mean_vector() -> np.ndarray | None:
    from phageforge.config import CACHE_DIR

    path = CACHE_DIR / MEAN_VECTOR_FILE
    return np.load(path) if path.exists() else None


def clean_sequence(protein: str) -> str:
    """Strip the trailing stop codon and normalise unknown residues."""
    seq = protein.rstrip("*").replace("*", "X")
    return seq.upper()


def to_documents(
    records: Iterable[ResolvedRBP],
    vectors: dict[str, np.ndarray],
    *,
    source: str,
    now: str,
) -> Iterator[tuple[str, dict]]:
    """Yield ``(protein_id, pf-proteins document)`` for every embedded RBP."""
    for record in records:
        if record.gene is None or record.protein_id not in vectors:
            continue
        sequence = clean_sequence(record.gene.protein)
        yield (
            record.protein_id,
            {
                "protein_id": record.protein_id,
                "locus_tag": record.rbp_id,
                "parent_id": record.phage_id,
                "parent_type": "phage",
                "role": "RBP",
                "rbp_type": record.rbp_type,
                "confidence": record.confidence,
                "stated_ordinal": record.stated_ordinal,
                "resolved_ordinal": record.resolved_ordinal,
                "ordinal_offset": record.offset,
                "resolution_note": record.note,
                "contig": record.gene.contig,
                "sequence": sequence,
                "length": len(sequence),
                "esm_vector": vectors[record.protein_id].tolist(),
                "source": source,
                "ingested_at": now,
            },
        )


def primary_rbp(records: Iterable[ResolvedRBP]) -> dict[str, ResolvedRBP]:
    """Pick one representative RBP per phage for the denormalised phage vector.

    Preference order: ``high`` over ``heuristic``, then the longest protein.
    Stage A needs a single vector per phage (ARCHITECTURE.md 7.1); which one it
    gets should be the best-evidenced, not an arbitrary first.
    """
    rank = {HIGH: 0, HEURISTIC: 1, UNRESOLVED: 2}
    best: dict[str, ResolvedRBP] = {}
    for record in records:
        if record.gene is None:
            continue
        current = best.get(record.phage_id)
        if current is None or (rank[record.confidence], -record.gene.length) < (
            rank[current.confidence],
            -current.gene.length,  # type: ignore[union-attr]
        ):
            best[record.phage_id] = record
    return best
