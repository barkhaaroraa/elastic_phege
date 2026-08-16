"""Index mappings for the six PhageForge indices.

Mappings are explicit and created up front. Dynamic mapping is disabled for the
vector fields in particular: a dynamically-mapped ``dense_vector`` is not indexed
for kNN, which fails silently at query time rather than at ingest.

Quantization (see ARCHITECTURE.md 5.4):
  * 256-dim genome sketches -> ``int8_hnsw`` (4x compression, low recall loss).
    BBQ needs higher dimensionality to be a good trade at this width.
  * protein embeddings      -> ``bbq_hnsw`` (~32x compression, shines at high dims).
"""

from __future__ import annotations

from elasticsearch import Elasticsearch

from phageforge.config import (
    BACTERIA,
    ESM_DIMS,
    INTERACTIONS,
    LOCAL_INDICES,
    PHAGES,
    PREDICTIONS,
    PROTEINS,
    SKETCH_DIMS,
)

# Single-node local dev: no replicas, or the cluster sits permanently yellow.
_LOCAL_SETTINGS = {"number_of_shards": 1, "number_of_replicas": 0}


#: Quantization for the 256-dim sketch vectors. **Deliberately off at this scale.**
#:
#: int8 quantization uses a single scale per vector. MDS coordinates decay sharply
#: in magnitude across dimensions, so the trailing dimensions quantize to zero and
#: the fine phylogenetic structure is destroyed -- measurably: with int8 a strain
#: did not retrieve *itself* as its own nearest neighbour, because quantization
#: error (~4e-3) exceeded typical inter-strain distance.
#:
#: At 404 strains x 256 dims, float32 costs ~413 KB. The saving is not worth the
#: accuracy. Quantization becomes the right call at the 100K-1M benchmark tiers,
#: where it is a headline variable to *measure* rather than assume.
SKETCH_INDEX_OPTIONS = {"type": "hnsw", "m": 16, "ef_construction": 100}


def _sketch_vector(similarity: str = "cosine") -> dict:
    """A 256-dim sketch vector.

    ``similarity`` differs by index, and the distinction is not cosmetic:

    * **phages** use tetranucleotide-frequency profiles -- compositional data,
      where cosine is the standard and correct measure.
    * **bacteria** use classical-MDS coordinates derived from the PanACoTA
      phylogenetic distance matrix. Those coordinates are *Euclidean by
      construction* -- the whole point of MDS is that Euclidean distance in the
      embedding reproduces the original distance. Scoring them by cosine would
      compare directions from an arbitrary origin and silently distort the
      phylogeny, so they use ``l2_norm``.
    """
    return {
        "type": "dense_vector",
        "dims": SKETCH_DIMS,
        "index": True,
        "similarity": similarity,
        "index_options": dict(SKETCH_INDEX_OPTIONS),
    }


def _protein_vector() -> dict:
    return {
        "type": "dense_vector",
        "dims": ESM_DIMS,
        "index": True,
        "similarity": "cosine",
        "index_options": {"type": "bbq_hnsw", "m": 16, "ef_construction": 100},
    }


MAPPINGS: dict[str, dict] = {
    PHAGES: {
        "settings": _LOCAL_SETTINGS,
        "mappings": {
            "properties": {
                "phage_id": {"type": "keyword"},
                "name": {"type": "text"},
                "accession": {"type": "keyword"},
                "taxonomy": {
                    "properties": {
                        "realm": {"type": "keyword"},
                        "family": {"type": "keyword"},
                        "subfamily": {"type": "keyword"},
                        "genus": {"type": "keyword"},
                        "species": {"type": "keyword"},
                    }
                },
                "morphotype": {"type": "keyword"},
                "lifestyle": {"type": "keyword"},
                "lifestyle_conf": {"type": "float"},
                "genome_length": {"type": "integer"},
                "host_range_summary": {"type": "text"},
                "known_host_species": {"type": "keyword"},
                "isolation_host": {"type": "keyword"},
                "isolation_host_phylo": {"type": "keyword"},
                "rbp_ids": {"type": "keyword"},
                "rbp_types": {"type": "keyword"},
                "has_rbp": {"type": "boolean"},
                # Which RBP the denormalised vector below came from, and how well
                # its locus tag resolved (see features/proteins.py).
                "rbp_match_id": {"type": "keyword"},
                "rbp_match_confidence": {"type": "keyword"},
                "receptor_classes": {"type": "keyword"},
                "has_depolymerase": {"type": "boolean"},
                "annotations": {"type": "text", "analyzer": "english"},
                "source": {"type": "keyword"},
                "ingested_at": {"type": "date"},
                # Breadth of known host range -- the "generalist phage" baseline
                # feature, precomputed so the ranker can learn it explicitly.
                "host_range_breadth": {"type": "float"},
                "genome_vector": _sketch_vector(),
                # Denormalised primary RBP embedding: avoids a protein->phage join
                # on every Stage A request.
                "rbp_match_vector": _protein_vector(),
            }
        },
    },
    BACTERIA: {
        "settings": _LOCAL_SETTINGS,
        "mappings": {
            "properties": {
                "strain_id": {"type": "keyword"},
                "species": {"type": "keyword"},
                "genus": {"type": "keyword"},
                "st": {"type": "keyword"},
                "phylogroup": {"type": "keyword"},
                # Surface/receptor biology. All keyword: Stage C runs
                # significant_terms over these, which needs doc values.
                "o_antigen": {"type": "keyword"},
                "h_type": {"type": "keyword"},
                "lps_type": {"type": "keyword"},
                "capsule_types": {"type": "keyword"},
                "abc_serotype": {"type": "keyword"},
                "k_locus": {"type": "keyword"},
                "defence_systems": {"type": "keyword"},
                "n_defense_systems": {"type": "integer"},
                "n_infections": {"type": "integer"},
                "pathotype": {"type": "keyword"},
                "isolation_source": {"type": "keyword"},
                "host_of_origin": {"type": "keyword"},
                "collection": {"type": "keyword"},
                "gembase": {"type": "keyword"},
                "source": {"type": "keyword"},
                "ingested_at": {"type": "date"},
                # Derived by MDS over the PanACoTA distance matrix -- the archive
                # ships no bacterial genomes to sketch. See DATA_NOTES.md.
                "genome_vector": _sketch_vector("l2_norm"),
            }
        },
    },
    INTERACTIONS: {
        "settings": _LOCAL_SETTINGS,
        "mappings": {
            "properties": {
                "phage_id": {"type": "keyword"},
                "host_id": {"type": "keyword"},
                "host_species": {"type": "keyword"},
                "infects": {"type": "boolean"},
                "eop": {"type": "float"},
                "score_raw": {"type": "float"},
                "assay_type": {"type": "keyword"},
                # strain | species. Load-bearing: only `strain` rows may be used
                # for training and evaluation.
                "resolution": {"type": "keyword"},
                "evidence_tier": {"type": "byte"},
                "source": {"type": "keyword"},
                "doi": {"type": "keyword"},
                "observed_at": {"type": "date"},
                "ingested_at": {"type": "date"},
            }
        },
    },
    PROTEINS: {
        "settings": _LOCAL_SETTINGS,
        "mappings": {
            "properties": {
                "protein_id": {"type": "keyword"},
                # The source locus tag. Not unique on its own -- one tag is
                # listed against two phages -- hence the composite protein_id.
                "locus_tag": {"type": "keyword"},
                "parent_id": {"type": "keyword"},
                "parent_type": {"type": "keyword"},
                "role": {"type": "keyword"},
                "rbp_type": {"type": "keyword"},
                # high | heuristic | unresolved -- how the locus tag was mapped
                # onto a called CDS. Recorded, never averaged away.
                "confidence": {"type": "keyword"},
                "stated_ordinal": {"type": "integer"},
                "resolved_ordinal": {"type": "integer"},
                "ordinal_offset": {"type": "integer"},
                # Why a non-`high` record ended up where it did. Provenance the
                # reviewer can read, not a number they have to trust.
                "resolution_note": {"type": "keyword"},
                "contig": {"type": "keyword"},
                "sequence": {"type": "text", "index": False},
                "length": {"type": "integer"},
                "source": {"type": "keyword"},
                "ingested_at": {"type": "date"},
                "esm_vector": _protein_vector(),
            }
        },
    },
    PREDICTIONS: {
        "settings": _LOCAL_SETTINGS,
        "mappings": {
            "properties": {
                "run_id": {"type": "keyword"},
                "strain_id": {"type": "keyword"},
                "created_at": {"type": "date"},
                "model_version": {"type": "keyword"},
                # Free-form per-stage detail (queries, counts, notes). Disabled
                # rather than dynamic: this is an audit blob, not a search surface.
                "stages": {"type": "object", "enabled": False},
                "shortlist": {"type": "object", "enabled": False},
                # Indexed summary fields, so the observability dashboard and the
                # /benchmark endpoint can aggregate runs without unpacking blobs.
                "shortlist_ids": {"type": "keyword"},
                "n_candidates": {"type": "integer"},
                "n_shortlist": {"type": "integer"},
                "stage_ms": {
                    "properties": {
                        "a": {"type": "float"},
                        "b": {"type": "float"},
                        "c": {"type": "float"},
                        "d": {"type": "float"},
                    }
                },
                "total_ms": {"type": "float"},
            }
        },
    },
}


def create_indices(
    client: Elasticsearch,
    indices: tuple[str, ...] = LOCAL_INDICES,
    *,
    recreate: bool = False,
    sync: bool = True,
) -> dict[str, str]:
    """Create the given indices. Returns ``{index: action}``.

    With ``recreate=True`` existing indices are deleted first -- destructive, so
    it is opt-in and never the default.

    With ``sync=True`` (the default) an index that already exists has its mapping
    PUT again. Elasticsearch adds new fields and rejects incompatible changes to
    existing ones, which is the behaviour we want: adding a field as the schema
    grows is routine and should not require a reindex, while a field whose type
    changed must fail loudly rather than be silently ignored.
    """
    actions: dict[str, str] = {}
    for name in indices:
        body = MAPPINGS[name]
        exists = client.indices.exists(index=name)
        if exists and recreate:
            client.indices.delete(index=name)
            exists = False
            actions[name] = "recreated"
        if exists:
            if sync:
                client.indices.put_mapping(index=name, **body["mappings"])
                actions.setdefault(name, "synced")
            else:
                actions.setdefault(name, "exists")
            continue
        client.indices.create(index=name, **body)
        actions.setdefault(name, "created")
    return actions
