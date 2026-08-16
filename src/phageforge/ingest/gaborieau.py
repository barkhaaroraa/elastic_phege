"""Parse the Gaborieau et al. (2024) archive into PhageForge documents.

See ``docs/DATA_NOTES.md`` for what the archive actually contains -- it differs
from the architecture doc in several places, and those notes are authoritative.

The three things worth knowing while reading this module:

* Every CSV in this archive is **semicolon-delimited**; the distance matrix and a
  few feature tables are tab-delimited.
* The interaction matrix holds a **graded 0-4 lysis score**, not a boolean and not
  EOP. We keep the raw score and binarise with a configurable threshold.
* Empty cells mean *untested*, which is not the same as *does not infect*. They
  are dropped, never coerced to negatives.
"""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from phageforge.ingest.zenodo import EXTRACT_DIR

ROOT = EXTRACT_DIR / "coli_phage_interactions_2023"
DATA = ROOT / "data"

MATRIX_CSV = DATA / "interactions" / "interaction_matrix.csv"
BACTERIA_CSV = DATA / "genomics" / "bacteria" / "picard_collection.csv"
DEFENCE_CSV = (
    DATA / "genomics" / "bacteria" / "defense_finder" / "370+host_defense_systems_subtypes.csv"
)
PHAGES_CSV = DATA / "genomics" / "phages" / "guelin_collection.csv"
RBP_CSV = DATA / "genomics" / "phages" / "RBP" / "RBP_list.csv"
FNA_DIR = DATA / "genomics" / "phages" / "FNA"
DISTANCE_TSV = (
    DATA
    / "genomics"
    / "bacteria"
    / "isolation_strains"
    / "panacota"
    / "tree"
    / "370+host_distance_matrix_orignames.tsv"
)

SOURCE = "gaborieau2024"
DOI = "10.1038/s41564-024-01832-5"

#: Lysis score above which we call an interaction positive. The authors' own
#: analysis code compares ``> 0`` throughout, so that is the default -- but this
#: is review question 3 in SOLUTION_OVERVIEW.md and is expected to move.
INFECT_THRESHOLD = float(os.environ.get("PF_INFECT_THRESHOLD", "0"))

_NULLS = {"", "NA", "na", "NaN", "nan", "-", "None"}


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _clean(value: str | None) -> str | None:
    """Normalise a raw CSV cell to a value or ``None``."""
    if value is None:
        return None
    value = value.strip()
    return None if value in _NULLS else value


def _read_delimited(path: Path, delimiter: str = ";") -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


# --------------------------------------------------------------- interactions


@dataclass
class MatrixStats:
    strains: int
    phages: int
    cells: int
    measured: int
    untested: int
    positives: int
    negatives: int
    score_counts: dict[str, int]


#: Host attributes copied from ``pf-bacteria`` onto every interaction row, as
#: ``bacteria field -> interaction field``. See the mapping comment in
#: ``es/mappings.py`` for why the denormalisation is worth its cost.
DENORMALISED_HOST_FIELDS = {
    "lps_type": "host_lps_type",
    "o_antigen": "host_o_antigen",
    "capsule_types": "host_capsule_types",
    "defence_systems": "host_defence_systems",
}


def read_interactions(hosts: dict[str, dict] | None = None) -> tuple[list[dict], MatrixStats]:
    """Parse the interaction matrix into flat interaction documents.

    ``hosts`` maps ``strain_id -> bacteria document`` and supplies the
    denormalised receptor/defence fields. It is read from :func:`read_bacteria`
    when not given; callers that already hold the bacteria records should pass
    them rather than paying for a second parse.
    """
    if hosts is None:
        hosts = {doc["strain_id"]: doc for doc in read_bacteria()}

    with MATRIX_CSV.open(newline="", encoding="utf-8-sig") as fh:
        rows = [line.rstrip("\n") for line in fh if line.strip()]

    header = rows[0].split(";")
    phage_ids = [p.strip() for p in header[1:]]

    docs: list[dict] = []
    score_counts: dict[str, int] = {}
    untested = 0
    now = _now()

    for row in rows[1:]:
        cells = row.split(";")
        strain_id = cells[0].strip()
        for phage_id, raw in zip(phage_ids, cells[1:], strict=False):
            raw = raw.strip()
            score_counts[raw or "(empty)"] = score_counts.get(raw or "(empty)", 0) + 1
            if raw in _NULLS:
                # Untested. Never coerce to a negative -- doing so would fabricate
                # ~157 measurements and quietly corrupt the benchmark.
                untested += 1
                continue
            score = float(raw)
            host = hosts.get(strain_id, {})
            denormalised = {
                target: host[source]
                for source, target in DENORMALISED_HOST_FIELDS.items()
                # A missing or null attribute is left off the document entirely
                # rather than written as null: Stage D reads these through terms
                # aggregations, where an absent field is correctly invisible and
                # an explicit null would still need special-casing.
                if host.get(source)
            }
            docs.append(
                {
                    "phage_id": phage_id,
                    "host_id": strain_id,
                    "host_species": "Escherichia coli",
                    "infects": score > INFECT_THRESHOLD,
                    **denormalised,
                    "score_raw": score,
                    "assay_type": "spot_lysis_graded",
                    "resolution": "strain",
                    "evidence_tier": 1,
                    "source": SOURCE,
                    "doi": DOI,
                    "ingested_at": now,
                }
            )

    positives = sum(1 for d in docs if d["infects"])
    stats = MatrixStats(
        strains=len(rows) - 1,
        phages=len(phage_ids),
        cells=(len(rows) - 1) * len(phage_ids),
        measured=len(docs),
        untested=untested,
        positives=positives,
        negatives=len(docs) - positives,
        score_counts=dict(sorted(score_counts.items())),
    )
    return docs, stats


# ------------------------------------------------------------------- bacteria


def _defence_systems() -> dict[str, list[str]]:
    """Map ``strain_id -> [defence system names present]``.

    The source table is wide: one column per system, cells are counts.
    """
    if not DEFENCE_CSV.exists():
        return {}
    out: dict[str, list[str]] = {}
    for row in _read_delimited(DEFENCE_CSV):
        strain = _clean(row.get("bacteria"))
        if not strain:
            continue
        systems = [
            name
            for name, count in row.items()
            if name != "bacteria" and _clean(count) not in (None, "0")
        ]
        out[strain] = sorted(systems)
    return out


_CAPSULE_COLUMNS = (
    "Capsule_ABC",
    "Capsule_GroupIV_e",
    "Capsule_GroupIV_e_stricte",
    "Capsule_GroupIV_s",
    "Capsule_Wzy_stricte",
)


def read_bacteria() -> list[dict]:
    """Parse strain metadata into ``pf-bacteria`` documents."""
    defences = _defence_systems()
    now = _now()
    docs: list[dict] = []

    for row in _read_delimited(BACTERIA_CSV):
        strain_id = _clean(row.get("bacteria"))
        if not strain_id:
            continue

        capsule_types = [
            col for col in _CAPSULE_COLUMNS if (_clean(row.get(col)) or "0") not in ("0", "0.0")
        ]

        docs.append(
            {
                "strain_id": strain_id,
                "species": "Escherichia coli",
                "genus": "Escherichia",
                "st": _clean(row.get("ST_Warwick")),
                "phylogroup": _clean(row.get("Clermont_Phylo")),
                "o_antigen": _clean(row.get("O-type")),
                "h_type": _clean(row.get("H-type")),
                "lps_type": _clean(row.get("LPS_type")),
                "capsule_types": capsule_types,
                "abc_serotype": _clean(row.get("ABC_serotype")),
                "pathotype": _clean(row.get("Pathotype")),
                "isolation_source": _clean(row.get("Origin")),
                "host_of_origin": _clean(row.get("Host")),
                "collection": _clean(row.get("Collection")),
                "gembase": _clean(row.get("Gembase")),
                "defence_systems": defences.get(strain_id, []),
                "n_defense_systems": _to_int(row.get("n_defense_systems")),
                "n_infections": _to_int(row.get("n_infections")),
                "source": SOURCE,
                "ingested_at": now,
            }
        )
    return docs


def _to_int(value: str | None) -> int | None:
    cleaned = _clean(value)
    if cleaned is None:
        return None
    try:
        return int(float(cleaned))
    except ValueError:
        return None


# --------------------------------------------------------------------- phages


def _rbp_index() -> dict[str, list[dict[str, str]]]:
    """Map ``phage_id -> [{rbp_id, type}]`` from ``RBP_list.csv``."""
    if not RBP_CSV.exists():
        return {}
    out: dict[str, list[dict[str, str]]] = {}
    for row in _read_delimited(RBP_CSV):
        phage = _clean(row.get("phage"))
        rbp_id = _clean(row.get("RBP"))
        if not phage or not rbp_id:
            continue
        out.setdefault(phage, []).append(
            {"rbp_id": rbp_id, "type": _clean(row.get("type")) or "unknown"}
        )
    return out


def _fna_lengths() -> dict[str, int]:
    """Genome length per phage, read from the FASTA files."""
    lengths: dict[str, int] = {}
    if not FNA_DIR.exists():
        return lengths
    for path in FNA_DIR.glob("*.fna"):
        total = 0
        with path.open() as fh:
            for line in fh:
                if not line.startswith(">"):
                    total += len(line.strip())
        lengths[path.stem] = total
    return lengths


def read_phages(interactions: list[dict] | None = None) -> list[dict]:
    """Parse phage metadata into ``pf-phages`` documents.

    ``host_range_breadth`` -- the fraction of tested strains a phage infects -- is
    computed here when interactions are supplied. It exists so the ranker can
    learn the "generalist phage" effect explicitly rather than being flattered by
    it (see TECHNICAL_DESIGN.md 7.4).
    """
    rbps = _rbp_index()
    lengths = _fna_lengths()
    now = _now()

    breadth: dict[str, float] = {}
    if interactions:
        tallies: dict[str, list[int]] = {}
        for rec in interactions:
            bucket = tallies.setdefault(rec["phage_id"], [0, 0])
            bucket[1] += 1
            if rec["infects"]:
                bucket[0] += 1
        breadth = {pid: hits / total for pid, (hits, total) in tallies.items() if total}

    docs: list[dict] = []
    for row in _read_delimited(PHAGES_CSV):
        phage_id = _clean(row.get("phage"))
        if not phage_id:
            continue
        entries = rbps.get(phage_id, [])
        declared = _to_int(row.get("Genome_size"))

        docs.append(
            {
                "phage_id": phage_id,
                "name": phage_id,
                "taxonomy": {
                    "family": _clean(row.get("Family")),
                    "genus": _clean(row.get("Genus")),
                    "subfamily": _clean(row.get("Subfamily")),
                    "species": _clean(row.get("Species")),
                },
                "morphotype": _clean(row.get("Morphotype")),
                # BACPHLIP is not run in this archive, so lifestyle is genuinely
                # unknown rather than assumed virulent.
                "lifestyle": "unknown",
                "genome_length": lengths.get(phage_id, declared),
                "isolation_host": _clean(row.get("Phage_host")),
                "isolation_host_phylo": _clean(row.get("Phage_host_phylo")),
                "known_host_species": ["Escherichia coli"],
                "rbp_ids": [e["rbp_id"] for e in entries],
                "rbp_types": sorted({e["type"] for e in entries}),
                "has_rbp": bool(entries),
                "host_range_breadth": breadth.get(phage_id),
                "annotations": " ".join(
                    filter(
                        None,
                        [
                            _clean(row.get("Family")),
                            _clean(row.get("Genus")),
                            _clean(row.get("Species")),
                            _clean(row.get("Morphotype")),
                            _clean(row.get("Subfamily")),
                            *[e["type"] for e in entries],
                        ],
                    )
                ),
                "source": SOURCE,
                "ingested_at": now,
            }
        )
    return docs


# ----------------------------------------------------------------- distances


def read_distance_matrix() -> tuple[list[str], list[list[float]]]:
    """Read the PanACoTA strain-strain distance matrix.

    Returns ``(strain_ids, matrix)``. Used to derive bacterial embeddings, since
    the archive ships no bacterial genomes to sketch (see DATA_NOTES.md).
    """
    with DISTANCE_TSV.open(newline="", encoding="utf-8-sig") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh if line.strip()]

    labels = [c.strip() for c in rows[0][1:]]
    ids: list[str] = []
    matrix: list[list[float]] = []
    for row in rows[1:]:
        ids.append(row[0].strip())
        matrix.append([float(v) if v.strip() not in _NULLS else 0.0 for v in row[1:]])

    if ids != labels:
        # Not fatal, but the matrix would be transposed relative to expectation.
        raise ValueError(
            f"distance matrix row/column labels disagree "
            f"({len(ids)} rows vs {len(labels)} cols; first mismatch matters)"
        )
    return ids, matrix
