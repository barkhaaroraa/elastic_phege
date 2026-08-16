"""The funnel stages. Every stage is an Elasticsearch operation.

Each function returns a :class:`StageResult` carrying results *and* metadata --
counts, timing, and the query that produced them. The metadata is what gets
persisted to ``pf-predictions`` and rendered in the UI, and it is the reason a
stage never simply returns an empty list: an empty stage raises
:class:`EmptyStage` naming what filtered everything out, because a blank
shortlist that looks like a confident "no phages" is the worst possible output.

**Stage ordering deviates from ARCHITECTURE.md 3**, which runs A (candidates)
before B (prior). With this dataset the dependency runs the other way. The
architecture's Stage A seeds its vector arms from the *strain's own* receptor
proteins, and the Gaborieau archive ships no bacterial protein sequences, so
there is nothing to embed. The adapted Stage A (DATA_NOTES.md 4) instead seeds
from the phages that already work on genomically nearby strains -- which is
exactly what Stage B computes. So the real order is:

    N  neighbours   kNN over pf-bacteria.genome_vector
    B  prior        weighted terms agg over pf-interactions   (seeds A)
    A  candidates   rrf: kNN(RBP) + kNN(genome) + BM25(annotations)
    C  explanation  significant_terms over pf-bacteria
    D  ranking      script_score blend over the survivors

At 96 phages Stage A cannot narrow anything -- the pool is larger than the
corpus. It is built and measured anyway because it is the stage that carries the
weight at the 100K-1M benchmark tiers, and because the ablation needs it.
"""

from __future__ import annotations

import math
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

from elasticsearch import Elasticsearch

from phageforge import config


class EmptyStage(RuntimeError):
    """A stage produced nothing, with a named cause.

    Raised rather than returning ``[]`` so the caller cannot mistake "the filter
    excluded everything" for "there is nothing to recommend".
    """

    def __init__(self, stage: str, reason: str) -> None:
        self.stage = stage
        self.reason = reason
        super().__init__(f"stage {stage} produced no results: {reason}")


@dataclass
class StageResult:
    name: str
    results: Any
    metadata: dict
    ms: float = 0.0


class _Timer:
    def __enter__(self):
        self.start = time.perf_counter()
        return self

    def __exit__(self, *exc):
        self.ms = (time.perf_counter() - self.start) * 1000


# ------------------------------------------------------------------ neighbours


@dataclass
class Neighbour:
    strain_id: str
    score: float
    distance: float
    weight: float
    fields: dict = field(default_factory=dict)


def _distance_from_score(score: float) -> float:
    """Invert Elasticsearch's ``l2_norm`` score back to Euclidean distance.

    ES scores ``l2_norm`` as ``1 / (1 + d^2)``. Even after the unit-RMS scaling
    applied to the MDS coordinates (see current_status.md 4.2), close strains
    still score 0.99996-1.0 -- the global range is usable but *within a
    neighbourhood* the scores are flat to four decimal places, so using them
    directly as similarity weights would weight all 25 neighbours equally and
    quietly turn Stage B into an unweighted average.

    Distance recovers the resolution exactly: the same neighbourhood spans
    d = 0.0 to 0.21. The kernel below then turns that into a real weight.
    """
    if score <= 0:
        return float("inf")
    return math.sqrt(max(1.0 / score - 1.0, 0.0))


def _kernel_weights(distances: Sequence[float]) -> list[float]:
    """Gaussian kernel with a bandwidth set from the neighbourhood itself.

    ``sigma`` is the median neighbour distance, so the weighting adapts to how
    dense the target's part of the phylogeny is: a strain in a tight ST cluster
    gets sharply peaked weights, an isolated one gets flatter ones. A fixed
    bandwidth would do the wrong thing at one end or the other.
    """
    finite = [d for d in distances if math.isfinite(d) and d > 0]
    if not finite:
        return [1.0] * len(distances)
    sigma = sorted(finite)[len(finite) // 2] or min(finite)
    return [math.exp(-((d / sigma) ** 2)) if math.isfinite(d) else 0.0 for d in distances]


NEIGHBOUR_FIELDS = (
    "strain_id",
    "st",
    "phylogroup",
    "o_antigen",
    "h_type",
    "lps_type",
    "capsule_types",
    "defence_systems",
)


def get_strain(client: Elasticsearch, strain_id: str) -> dict:
    """Fetch a strain with its vector. Raises if it is missing or unvectorised."""
    hits = client.search(
        index=config.BACTERIA,
        size=1,
        query={"term": {"strain_id": strain_id}},
        source_includes=[*NEIGHBOUR_FIELDS, "genome_vector", "species"],
    )["hits"]["hits"]
    if not hits:
        raise EmptyStage("N", f"strain {strain_id!r} is not in {config.BACTERIA}")
    source = hits[0]["_source"]
    if not source.get("genome_vector"):
        # LF110 is the known case: indexed, but absent from the distance matrix.
        raise EmptyStage("N", f"strain {strain_id!r} has no genome_vector")
    return source


def find_neighbours(
    client: Elasticsearch,
    strain: dict,
    *,
    k: int = config.NEIGHBOUR_K,
    exclude: Sequence[str] = (),
) -> StageResult:
    """kNN over ``pf-bacteria.genome_vector`` for the k most similar strains.

    ``exclude`` must contain the target strain itself. It is *not* excluded
    implicitly, because the benchmark also needs to exclude a whole held-out
    fold, and a stage that silently drops one document is harder to reason about
    than one that is told what to drop.
    """
    target = strain["strain_id"]
    excluded = {target, *exclude}
    with _Timer() as t:
        response = client.search(
            index=config.BACTERIA,
            size=k + len(excluded),
            knn={
                "field": "genome_vector",
                "query_vector": strain["genome_vector"],
                "k": k + len(excluded),
                "num_candidates": max(500, (k + len(excluded)) * 10),
                "filter": {"term": {"species": strain.get("species", "Escherichia coli")}},
            },
            source_includes=list(NEIGHBOUR_FIELDS),
        )

    hits = [h for h in response["hits"]["hits"] if h["_source"]["strain_id"] not in excluded][:k]
    if not hits:
        raise EmptyStage("N", f"no tested strains within kNN of {target!r} after exclusions")

    distances = [_distance_from_score(h["_score"]) for h in hits]
    weights = _kernel_weights(distances)
    neighbours = [
        Neighbour(
            strain_id=h["_source"]["strain_id"],
            score=float(h["_score"]),
            distance=d,
            weight=w,
            fields=h["_source"],
        )
        for h, d, w in zip(hits, distances, weights, strict=True)
    ]
    return StageResult(
        name="N",
        results=neighbours,
        metadata={
            "k_requested": k,
            "n_returned": len(neighbours),
            "excluded": sorted(excluded),
            "distance_min": min(distances),
            "distance_max": max(distances),
            "weight_min": min(weights),
            "weight_max": max(weights),
            "same_st": sum(
                1 for n in neighbours if n.fields.get("st") and n.fields["st"] == strain.get("st")
            ),
        },
        ms=t.ms,
    )


# ---------------------------------------------------------------- B: the prior


@dataclass
class Prior:
    phage_id: str
    prior: float  # similarity-weighted infection rate across neighbours
    hit_rate: float  # unweighted, for comparison
    n_tested: int
    n_hits: int
    weight_sum: float


def stage_b_prior(
    client: Elasticsearch,
    neighbours: Sequence[Neighbour],
    *,
    size: int = 2_000,
) -> StageResult:
    """Susceptibility transfer from genomically near strains.

    Collaborative filtering done entirely in Elasticsearch: one ``terms``
    aggregation over ``pf-interactions``, with a runtime field injecting each
    neighbour's kernel weight so the ``weighted_avg`` happens in-cluster rather
    than in pandas (ARCHITECTURE.md 7.2, TECHNICAL_DESIGN.md R3).

    Only ``resolution: strain`` rows count. Species-level evidence would inject
    label noise at exactly the resolution we claim to operate at.
    """
    weights = {n.strain_id: n.weight for n in neighbours}
    with _Timer() as t:
        response = client.search(
            index=config.INTERACTIONS,
            size=0,
            runtime_mappings={
                "neighbour_similarity": {
                    "type": "double",
                    "script": {
                        "source": (
                            "String h = doc['host_id'].value; "
                            "emit(params.w.containsKey(h) ? params.w.get(h) : 0.0)"
                        ),
                        "params": {"w": weights},
                    },
                }
            },
            query={
                "bool": {
                    "filter": [
                        {"terms": {"host_id": list(weights)}},
                        {"term": {"resolution": "strain"}},
                    ]
                }
            },
            aggs={
                "by_phage": {
                    "terms": {"field": "phage_id", "size": size},
                    "aggs": {
                        "hit_rate": {"avg": {"field": "infects"}},
                        "n_hits": {"sum": {"field": "infects"}},
                        "n_tested": {"value_count": {"field": "host_id"}},
                        "weight_sum": {"sum": {"field": "neighbour_similarity"}},
                        "weighted": {
                            "weighted_avg": {
                                "value": {"field": "infects"},
                                "weight": {"field": "neighbour_similarity"},
                            }
                        },
                    },
                }
            },
        )

    buckets = response["aggregations"]["by_phage"]["buckets"]
    if not buckets:
        raise EmptyStage(
            "B",
            f"none of the {len(weights)} neighbour strains has a strain-resolution "
            "interaction row",
        )

    priors = [
        Prior(
            phage_id=b["key"],
            prior=float(b["weighted"]["value"] or 0.0),
            hit_rate=float(b["hit_rate"]["value"] or 0.0),
            n_tested=int(b["n_tested"]["value"]),
            n_hits=int(b["n_hits"]["value"]),
            weight_sum=float(b["weight_sum"]["value"] or 0.0),
        )
        for b in buckets
    ]
    priors.sort(key=lambda p: -p.prior)
    return StageResult(
        name="B",
        results=priors,
        metadata={
            "n_phages": len(priors),
            "n_neighbours": len(weights),
            "n_with_positive_prior": sum(1 for p in priors if p.prior > 0),
            "top": [(p.phage_id, round(p.prior, 3)) for p in priors[:5]],
        },
        ms=t.ms,
    )


# ------------------------------------------------------------ A: candidate set


#: How many of Stage B's top phages seed Stage A's query vectors. These are the
#: phages that demonstrably work in this part of the phylogeny; Stage A's job is
#: to find phages *like* them that have never been tested here.
SEED_SIZE = 10


def _weighted_centroid(vectors: list[list[float]], weights: list[float]) -> list[float] | None:
    """L2-normalised weighted mean. ``None`` if there is nothing to average."""
    if not vectors:
        return None
    total = sum(weights) or 1.0
    dims = len(vectors[0])
    acc = [0.0] * dims
    for vec, w in zip(vectors, weights, strict=True):
        for i, v in enumerate(vec):
            acc[i] += v * w / total
    norm = math.sqrt(sum(v * v for v in acc))
    if norm == 0:
        return None
    return [v / norm for v in acc]


def _knn_scores(
    client: Elasticsearch, field_name: str, vector: list[float], k: int
) -> dict[str, float]:
    """Raw kNN similarity per phage for one arm.

    ES scores ``cosine`` as ``(1 + cos) / 2``, so this lands in [0, 1] with 0.5
    meaning orthogonal -- usable directly as a feature without rescaling.
    """
    response = client.search(
        index=config.PHAGES,
        size=k,
        knn={
            "field": field_name,
            "query_vector": vector,
            "k": k,
            "num_candidates": max(k * 2, 100),
        },
        source_includes=["phage_id"],
    )
    return {h["_source"]["phage_id"]: float(h["_score"]) for h in response["hits"]["hits"]}


def stage_a_candidates(
    client: Elasticsearch,
    priors: Sequence[Prior],
    *,
    pool: int = config.CANDIDATE_POOL,
    seed_size: int = SEED_SIZE,
    use_rbp: bool = True,
) -> StageResult:
    """Candidate generation: one ``rrf`` retriever fusing three signals.

    The three arms are deliberately independent, and RRF fuses them by rank so no
    score normalisation is needed across a BM25 score and two cosine similarities:

    1. **RBP similarity** -- kNN over ``rbp_match_vector`` against the ESM-2
       centroid of the seed phages' RBPs. This is the biologically load-bearing
       arm: adsorption dominates, so a phage whose fibre resembles a fibre that
       already works here is a real lead.
    2. **Genome composition** -- kNN over ``genome_vector`` (TNF) against the seed
       centroid. Cheap, and it recovers taxonomy (96.8% genus agreement).
    3. **Annotations** -- BM25 over the seed phages' taxonomy and RBP-type terms.
       A lexical arm catches relatives the vectors miss.

    ``use_rbp=False`` drops arm 1 -- that is the "-- RBP features" ablation in
    TECHNICAL_DESIGN.md 7.6, not a fallback.
    """
    seeds = [p for p in priors if p.prior > 0][:seed_size]
    if not seeds:
        raise EmptyStage(
            "A",
            "no phage infects any neighbour strain, so there is nothing to seed "
            "candidate generation with",
        )
    seed_weights = [p.prior for p in seeds]
    seed_ids = [p.phage_id for p in seeds]

    with _Timer() as t:
        seed_docs = client.search(
            index=config.PHAGES,
            size=len(seed_ids),
            query={"terms": {"phage_id": seed_ids}},
            source_includes=["phage_id", "genome_vector", "rbp_match_vector", "annotations"],
        )["hits"]["hits"]
        by_id = {h["_source"]["phage_id"]: h["_source"] for h in seed_docs}

        genome_vecs, genome_ws = [], []
        rbp_vecs, rbp_ws = [], []
        terms: list[str] = []
        for pid, w in zip(seed_ids, seed_weights, strict=True):
            source = by_id.get(pid)
            if not source:
                continue
            if source.get("genome_vector"):
                genome_vecs.append(source["genome_vector"])
                genome_ws.append(w)
            if use_rbp and source.get("rbp_match_vector"):
                rbp_vecs.append(source["rbp_match_vector"])
                rbp_ws.append(w)
            if source.get("annotations"):
                terms.append(source["annotations"])

        arms: list[dict] = []
        arm_names: list[str] = []
        rbp_centroid = _weighted_centroid(rbp_vecs, rbp_ws)
        if rbp_centroid:
            arms.append(
                {
                    "knn": {
                        "field": "rbp_match_vector",
                        "query_vector": rbp_centroid,
                        "k": pool,
                        "num_candidates": max(pool * 2, 100),
                    }
                }
            )
            arm_names.append("rbp")
        genome_centroid = _weighted_centroid(genome_vecs, genome_ws)
        if genome_centroid:
            arms.append(
                {
                    "knn": {
                        "field": "genome_vector",
                        "query_vector": genome_centroid,
                        "k": pool,
                        "num_candidates": max(pool * 2, 100),
                    }
                }
            )
            arm_names.append("genome")
        if terms:
            arms.append({"standard": {"query": {"match": {"annotations": " ".join(terms)}}}})
            arm_names.append("bm25")

        if not arms:
            raise EmptyStage("A", "seed phages carry neither vectors nor annotations")

        if len(arms) == 1:
            body: dict = {"retriever": arms[0]}
        else:
            body = {
                "retriever": {
                    "rrf": {
                        "retrievers": arms,
                        "rank_window_size": pool,
                        "rank_constant": 60,
                    }
                }
            }

        response = client.search(
            index=config.PHAGES,
            size=pool,
            source_includes=["phage_id", "host_range_breadth", "rbp_match_confidence"],
            **body,
        )

        # RRF fuses by *rank*, which is what makes it robust across heterogeneous
        # retrievers -- but it also discards the per-arm similarity, and Stage D
        # needs those as features or the RBP signal never reaches the score at
        # all (measurably: without this the "- RBP arm" ablation moved P@10 by
        # 0.000). So each vector arm is also run on its own to recover its score.
        arm_scores = {
            name: _knn_scores(client, field_name, vector, pool)
            for name, field_name, vector in (
                ("rbp", "rbp_match_vector", rbp_centroid),
                ("genome", "genome_vector", genome_centroid),
            )
            if vector is not None
        }

    hits = response["hits"]["hits"]
    if not hits:
        raise EmptyStage("A", "the rrf retriever matched no phages")
    candidates = [
        {
            "phage_id": h["_source"]["phage_id"],
            "rank_score": float(h["_score"]),
            "host_range_breadth": h["_source"].get("host_range_breadth"),
            "rbp_confidence": h["_source"].get("rbp_match_confidence"),
            "rbp_similarity": arm_scores.get("rbp", {}).get(h["_source"]["phage_id"], 0.0),
            "genome_similarity": arm_scores.get("genome", {}).get(h["_source"]["phage_id"], 0.0),
        }
        for h in hits
    ]
    return StageResult(
        name="A",
        results=candidates,
        metadata={
            "arms": arm_names,
            "seed_phages": seed_ids,
            "pool_requested": pool,
            "n_candidates": len(candidates),
            "corpus_size": int(client.count(index=config.PHAGES)["count"]),
        },
        ms=t.ms,
    )


# ------------------------------------------------------------- C: explanations


#: Bacterial features Stage C tests for enrichment. All keyword-mapped, because
#: significant_terms needs doc values.
EXPLAIN_FIELDS = (
    "o_antigen",
    "lps_type",
    "capsule_types",
    "st",
    "phylogroup",
    "defence_systems",
    "h_type",
)

_FIELD_LABELS = {
    "o_antigen": "O-antigen",
    "lps_type": "LPS core type",
    "capsule_types": "capsule",
    "st": "sequence type ST",
    "phylogroup": "phylogroup",
    "defence_systems": "defence system",
    "h_type": "H-type",
}


def stage_c_explain(
    client: Elasticsearch,
    phage_ids: Sequence[str],
    strain: dict,
    *,
    per_field: int = 3,
    max_susceptible: int = 500,
) -> StageResult:
    """Why would this phage infect this strain? ``significant_terms``, per phage.

    Among the strains a phage is known to infect, which features are *enriched*
    against the whole population -- not merely common. A feature present in 90%
    of susceptible strains and 90% of everything explains nothing;
    ``significant_terms`` scores exactly that difference, in-cluster.

    Reasons are then filtered to features the *target strain actually has*, which
    is what turns an enrichment table into "your strain is O25, and this phage
    infects O25 strains far above background".
    """
    if not phage_ids:
        raise EmptyStage("C", "no phages to explain")

    with _Timer() as t:
        # One request for the susceptible-strain sets of every phage at once.
        susceptible = client.search(
            index=config.INTERACTIONS,
            size=0,
            query={
                "bool": {
                    "filter": [
                        {"terms": {"phage_id": list(phage_ids)}},
                        {"term": {"resolution": "strain"}},
                        {"term": {"infects": True}},
                    ]
                }
            },
            aggs={
                "by_phage": {
                    "terms": {"field": "phage_id", "size": len(phage_ids)},
                    "aggs": {"hosts": {"terms": {"field": "host_id", "size": max_susceptible}}},
                }
            },
        )["aggregations"]["by_phage"]["buckets"]

        hosts_by_phage = {
            b["key"]: [h["key"] for h in b["hosts"]["buckets"]] for b in susceptible
        }

        # One msearch for every phage's significant_terms profile.
        searches: list[dict] = []
        ordered: list[str] = []
        for phage_id in phage_ids:
            hosts = hosts_by_phage.get(phage_id)
            if not hosts:
                continue
            ordered.append(phage_id)
            searches.append({"index": config.BACTERIA})
            searches.append(
                {
                    "size": 0,
                    "query": {"terms": {"strain_id": hosts}},
                    "aggs": {
                        f"sig_{field}": {
                            "significant_terms": {"field": field, "size": per_field}
                        }
                        for field in EXPLAIN_FIELDS
                    },
                }
            )
        profiles = client.msearch(searches=searches)["responses"] if searches else []

    explanations: dict[str, dict] = {}
    for phage_id, response in zip(ordered, profiles, strict=True):
        if "error" in response:
            explanations[phage_id] = {"reasons": [], "error": str(response["error"])[:200]}
            continue
        signals = []
        foreground = response["hits"]["total"]["value"]
        for field_name in EXPLAIN_FIELDS:
            agg = response["aggregations"][f"sig_{field_name}"]
            # Background size is reported once per aggregation, not per bucket.
            background = agg.get("bg_count") or 0
            for bucket in agg["buckets"]:
                target_value = strain.get(field_name)
                matches = (
                    bucket["key"] in target_value
                    if isinstance(target_value, list)
                    else bucket["key"] == target_value
                )
                bg_frac = bucket["bg_count"] / background if background else 0.0
                signals.append(
                    {
                        "field": field_name,
                        "value": bucket["key"],
                        "score": round(float(bucket["score"]), 4),
                        "fg_count": bucket["doc_count"],
                        "fg_total": foreground,
                        "bg_count": bucket["bg_count"],
                        "bg_total": background,
                        "fg_frac": bucket["doc_count"] / max(foreground, 1),
                        "bg_frac": bg_frac,
                        "matches_target": matches,
                    }
                )
        signals.sort(key=lambda s: (-int(s["matches_target"]), -s["score"]))
        explanations[phage_id] = {
            "n_susceptible": len(hosts_by_phage.get(phage_id, [])),
            "signals": signals,
            "reasons": [_reason(s) for s in signals if s["matches_target"]][:3],
        }

    n_with_reason = sum(1 for e in explanations.values() if e.get("reasons"))
    return StageResult(
        name="C",
        results=explanations,
        metadata={
            "n_phages": len(explanations),
            "n_with_matching_reason": n_with_reason,
            "fields": list(EXPLAIN_FIELDS),
        },
        ms=t.ms,
    )


#: Fields where the strain *carries* the feature rather than *is* it. Purely a
#: wording concern, but "your strain is PsyrTA" reads as nonsense to a
#: microbiologist and the reason strings are the part a reviewer judges.
_CARRIED_FIELDS = {"defence_systems", "capsule_types"}


def _reason(signal: dict) -> str:
    label = _FIELD_LABELS.get(signal["field"], signal["field"])
    verb = "carries" if signal["field"] in _CARRIED_FIELDS else "is"
    return (
        f"{signal['fg_count']}/{signal['fg_total']} of its known hosts are {label} "
        f"{signal['value']} vs {signal['bg_frac']:.0%} of all strains "
        f"(enrichment {signal['score']:.2f}); your strain {verb} {signal['value']}"
    )


# ----------------------------------------------------------------- D: ranking


#: Feature weights for the Stage D blend. These are the *untrained* defaults --
#: the neighbour prior dominates because it is the only feature measured against
#: real negatives. `phageforge.bench` fits them on strain-grouped folds and
#: writes the fitted set to data/derived; see MODEL_VERSION below.
DEFAULT_WEIGHTS = {
    "prior": 3.0,
    "support": 0.30,
    "candidate_rank": 0.50,
    "breadth": 0.25,
    "rbp_similarity": 0.0,
    "genome_similarity": 0.0,
    "bias": 0.0,
}

#: Every feature the Stage D script consumes, in a fixed order. The trainer fits
#: exactly these, so adding a feature here without retraining leaves it at its
#: default weight rather than silently changing the model.
FEATURE_NAMES = (
    "prior",
    "support",
    "candidate_rank",
    "breadth",
    "rbp_similarity",
    "genome_similarity",
)

MODEL_VERSION = "phageforge-linear-v1"


def stage_d_rank(
    client: Elasticsearch,
    candidates: Sequence[dict],
    priors: Sequence[Prior],
    *,
    top_n: int = config.TOP_N,
    weights: dict[str, float] | None = None,
    breadth: dict[str, float] | None = None,
) -> StageResult:
    """Final ranking: a linear blend evaluated in-cluster by ``script_score``.

    ARCHITECTURE.md 7.4 specifies a ``learning_to_rank`` rescorer running an
    XGBRanker pushed in via Eland. That needs a trained model, which needs the
    Phase 6 benchmark harness, so this is the same shape one step earlier: the
    features are computed by the previous stages and passed as script params,
    and the *scoring* runs inside Elasticsearch. Swapping this for the LTR
    rescorer is then a change of one query body, not of the architecture.

    Keeping the arithmetic in ES rather than sorting a Python list is the point:
    TECHNICAL_DESIGN.md R3 makes "Elastic becomes a metadata store" a fatal
    failure mode, and Stage D is where that would happen first.

    ``breadth`` overrides the indexed ``host_range_breadth``. The benchmark needs
    this: the indexed value is computed over *all* interactions including the
    held-out strain's own row, which is a small but real leak, so each fold
    supplies breadth computed from its training strains only.
    """
    if not candidates:
        raise EmptyStage("D", "stage A returned no candidates to rank")
    weights = {**DEFAULT_WEIGHTS, **(weights or {})}

    prior_by_id = {p.phage_id: p for p in priors}
    n = len(candidates)
    features = {
        c["phage_id"]: {
            "prior": prior_by_id[c["phage_id"]].prior if c["phage_id"] in prior_by_id else 0.0,
            # log1p of neighbour support: the difference between 1 and 5 tested
            # neighbours matters far more than between 20 and 24.
            "support": math.log1p(
                prior_by_id[c["phage_id"]].n_tested if c["phage_id"] in prior_by_id else 0
            ),
            # Stage A's fused rank, normalised so arm scores of different scales
            # cannot dominate; 1.0 for the best candidate, ~0 for the last.
            "candidate_rank": 1.0 - (i / n),
            "breadth": (
                float(breadth.get(c["phage_id"], 0.0))
                if breadth is not None
                else float(c.get("host_range_breadth") or 0.0)
            ),
        }
        for i, c in enumerate(candidates)
    }

    script = """
        String id = doc['phage_id'].value;
        if (!params.f.containsKey(id)) { return 0; }
        Map f = params.f.get(id);
        double breadth = (double) f.breadth;
        double s = params.w.bias
            + params.w.prior          * ((double) f.prior)
            + params.w.support        * ((double) f.support)
            + params.w.candidate_rank * ((double) f.candidate_rank)
            + params.w.breadth        * breadth;
        // script_score requires a non-negative score.
        return s < 0 ? 0 : s;
    """

    with _Timer() as t:
        response = client.search(
            index=config.PHAGES,
            size=top_n,
            query={
                "script_score": {
                    "query": {"terms": {"phage_id": list(features)}},
                    "script": {"source": script, "params": {"f": features, "w": weights}},
                }
            },
            source_includes=[
                "phage_id",
                "taxonomy",
                "morphotype",
                "lifestyle",
                "genome_length",
                "host_range_breadth",
                "rbp_ids",
                "rbp_match_confidence",
            ],
        )

    hits = response["hits"]["hits"]
    if not hits:
        raise EmptyStage("D", "script_score matched none of the candidate phage ids")

    ranked = []
    for h in hits:
        source = h["_source"]
        pid = source["phage_id"]
        prior = prior_by_id.get(pid)
        ranked.append(
            {
                "phage_id": pid,
                "score": float(h["_score"]),
                "taxonomy": source.get("taxonomy", {}),
                "morphotype": source.get("morphotype"),
                "lifestyle": source.get("lifestyle"),
                "genome_length": source.get("genome_length"),
                "host_range_breadth": source.get("host_range_breadth"),
                "rbp_confidence": source.get("rbp_match_confidence"),
                # Tier 1 whenever the recommendation rests on measured
                # strain-level interactions among the neighbours; Tier 4 when it
                # rests on vector similarity alone.
                "evidence_tier": 1 if prior and prior.n_tested else 4,
                "prior": prior.prior if prior else 0.0,
                "n_neighbours_tested": prior.n_tested if prior else 0,
                "n_neighbours_infected": prior.n_hits if prior else 0,
                "features": features[pid],
            }
        )
    return StageResult(
        name="D",
        results=ranked,
        metadata={
            "model_version": MODEL_VERSION,
            "weights": weights,
            "n_scored": len(features),
            "n_returned": len(ranked),
        },
        ms=t.ms,
    )
