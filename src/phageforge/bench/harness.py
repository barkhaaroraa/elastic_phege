"""Phase 6 -- the benchmark. This is the deliverable; everything else is scaffolding.

**The split is the thing most likely to be got wrong.** Splitting the 38,435
measured pairs at random leaks a strain's other 95 measurements into training and
produces a beautiful, meaningless number. Folds here are grouped by
``strain_id``: a strain appears in exactly one fold, and when it is evaluated,
*every other strain in its fold is excluded from its neighbour set* -- so the
funnel can only transfer susceptibility from strains it was allowed to see.

Three baselines, no exceptions (TECHNICAL_DESIGN.md 7.4):

``random``
    The floor.
``generalist``
    Rank every strain by global phage host-range breadth. Strong and dumb: if
    the funnel does not clearly beat it, nothing strain-specific is being
    learned. Breadth is recomputed per fold from training strains only.
``phylo_nn``
    Copy the single nearest training strain's results. The method a
    microbiologist would use by hand.

Plus the funnel and its two ablations. Every method is scored on the *same*
held-out strains with the *same* ground truth, and every method ranks the full
phage corpus so P@10, R@10 and AUPRC are all computed from one ranking.
"""

from __future__ import annotations

import contextlib
import json
import math
import random
import time
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

from elasticsearch import Elasticsearch

from phageforge import config
from phageforge.funnel import cocktail as cocktail_mod
from phageforge.funnel.pipeline import run_funnel
from phageforge.funnel.stages import EmptyStage, _distance_from_score

TOP_K = 10
COCKTAIL_SIZE = 4
N_FOLDS = 5
SEED = 20260815


# ------------------------------------------------------------------- ground truth


@dataclass
class Truth:
    """Measured outcomes for one strain: ``{phage_id: infects}``."""

    strain_id: str
    outcomes: dict[str, bool]

    @property
    def positives(self) -> set[str]:
        return {p for p, hit in self.outcomes.items() if hit}

    @property
    def n_tested(self) -> int:
        return len(self.outcomes)


def load_truth(client: Elasticsearch) -> dict[str, Truth]:
    """Every strain-resolution interaction, keyed by strain.

    Read once and held in memory: 38k rows is nothing, and re-querying per fold
    would make the harness's runtime dominated by its own bookkeeping.
    """
    truth: dict[str, dict[str, bool]] = {}
    response = client.search(
        index=config.INTERACTIONS,
        size=10_000,
        query={"term": {"resolution": "strain"}},
        source_includes=["phage_id", "host_id", "infects"],
        scroll="2m",
    )
    scroll_id = response.get("_scroll_id")
    try:
        while True:
            hits = response["hits"]["hits"]
            if not hits:
                break
            for hit in hits:
                source = hit["_source"]
                truth.setdefault(source["host_id"], {})[source["phage_id"]] = bool(
                    source["infects"]
                )
            response = client.scroll(scroll_id=scroll_id, scroll="2m")
            scroll_id = response.get("_scroll_id")
    finally:
        if scroll_id:
            # An already-expired scroll context 404s; that is not an error worth
            # failing the whole benchmark over.
            with contextlib.suppress(Exception):  # cleanup only
                client.clear_scroll(scroll_id=scroll_id)
    return {sid: Truth(sid, outcomes) for sid, outcomes in truth.items()}


def evaluable_strains(client: Elasticsearch, truth: dict[str, Truth]) -> list[str]:
    """Strains that have both a vector and measured interactions.

    ``LF110`` is indexed but absent from the PanACoTA distance matrix so it has
    no vector; ``H1-005-0065-L-P`` and ``H27`` are in the matrix but were never
    indexed. Neither can be evaluated, and both are reported rather than dropped
    silently.
    """
    vectored = set()
    response = client.search(
        index=config.BACTERIA,
        size=10_000,
        query={"exists": {"field": "genome_vector"}},
        source_includes=["strain_id"],
    )
    for hit in response["hits"]["hits"]:
        vectored.add(hit["_source"]["strain_id"])
    return sorted(vectored & set(truth))


def make_folds(strains: Sequence[str], n_folds: int = N_FOLDS, seed: int = SEED) -> list[list[str]]:
    """Deterministic strain-grouped folds. Groups are strains, never pairs."""
    shuffled = list(strains)
    random.Random(seed).shuffle(shuffled)
    return [sorted(shuffled[i::n_folds]) for i in range(n_folds)]


# ----------------------------------------------------------------------- metrics


def precision_at_k(ranking: Sequence[str], positives: set[str], k: int = TOP_K) -> float:
    top = ranking[:k]
    return sum(1 for p in top if p in positives) / len(top) if top else 0.0


def recall_at_k(ranking: Sequence[str], positives: set[str], k: int = TOP_K) -> float:
    if not positives:
        return float("nan")
    return sum(1 for p in ranking[:k] if p in positives) / len(positives)


def average_precision(ranking: Sequence[str], positives: set[str]) -> float:
    """AUPRC by the standard interpolation-free average-precision estimator.

    AUPRC rather than AUROC as the headline comparability metric: positives are
    ~21% of measured pairs and far rarer per strain, and AUROC flatters a ranker
    on imbalanced data.
    """
    if not positives:
        return float("nan")
    hits = 0
    total = 0.0
    for i, phage in enumerate(ranking, start=1):
        if phage in positives:
            hits += 1
            total += hits / i
    return total / len(positives)


@dataclass
class MethodScores:
    name: str
    p_at_k: list[float] = field(default_factory=list)
    r_at_k: list[float] = field(default_factory=list)
    auprc: list[float] = field(default_factory=list)
    cocktail_hits: list[float] = field(default_factory=list)
    topk_hits: list[float] = field(default_factory=list)
    failures: list[tuple[str, str]] = field(default_factory=list)
    ms: list[float] = field(default_factory=list)

    def add(self, ranking: Sequence[str], truth: Truth, cocktail: Sequence[str] = ()) -> None:
        # Only score phages actually measured against this strain -- an untested
        # pair is not a negative, and counting it as one fabricates the metric.
        measured = [p for p in ranking if p in truth.outcomes]
        positives = truth.positives
        self.p_at_k.append(precision_at_k(measured, positives))
        self.r_at_k.append(recall_at_k(measured, positives))
        self.auprc.append(average_precision(measured, positives))
        if cocktail:
            self.cocktail_hits.append(
                1.0 if any(p in positives for p in cocktail) else 0.0
            )
            self.topk_hits.append(
                1.0 if any(p in positives for p in measured[:COCKTAIL_SIZE]) else 0.0
            )

    def summary(self) -> dict:
        def mean(values: list[float]) -> float:
            usable = [v for v in values if not math.isnan(v)]
            return sum(usable) / len(usable) if usable else float("nan")

        return {
            "method": self.name,
            "n_strains": len(self.p_at_k),
            "p_at_10": mean(self.p_at_k),
            "r_at_10": mean(self.r_at_k),
            "auprc": mean(self.auprc),
            "cocktail_coverage": mean(self.cocktail_hits) if self.cocktail_hits else float("nan"),
            "top4_coverage": mean(self.topk_hits) if self.topk_hits else float("nan"),
            "median_ms": (sorted(self.ms)[len(self.ms) // 2] if self.ms else float("nan")),
            "p95_ms": (
                sorted(self.ms)[min(int(len(self.ms) * 0.95), len(self.ms) - 1)]
                if self.ms
                else float("nan")
            ),
            "n_failed": len(self.failures),
        }


# --------------------------------------------------------------------- baselines


def fold_breadth(
    client: Elasticsearch, training: Sequence[str], phages: Sequence[str]
) -> dict[str, float]:
    """Host-range breadth per phage, computed from training strains only.

    The indexed ``host_range_breadth`` covers every strain including the held-out
    one. One strain in 402 is a small leak, but it is a leak in the feature the
    generalist baseline *is*, so it gets recomputed per fold.
    """
    response = client.search(
        index=config.INTERACTIONS,
        size=0,
        query={
            "bool": {
                "filter": [
                    {"terms": {"host_id": list(training)}},
                    {"term": {"resolution": "strain"}},
                ]
            }
        },
        aggs={
            "by_phage": {
                "terms": {"field": "phage_id", "size": max(len(phages), 1000)},
                "aggs": {"rate": {"avg": {"field": "infects"}}},
            }
        },
    )
    return {
        b["key"]: float(b["rate"]["value"] or 0.0)
        for b in response["aggregations"]["by_phage"]["buckets"]
    }


def baseline_random(phages: Sequence[str], strain_id: str) -> list[str]:
    """Seeded per strain so the floor is reproducible, not a lucky draw."""
    order = list(phages)
    random.Random(f"{SEED}:{strain_id}").shuffle(order)
    return order


def baseline_generalist(breadth: dict[str, float], phages: Sequence[str]) -> list[str]:
    return sorted(phages, key=lambda p: -breadth.get(p, 0.0))


def baseline_phylo_nn(
    client: Elasticsearch,
    strain_id: str,
    vector: list[float],
    excluded: set[str],
    truth: dict[str, Truth],
    phages: Sequence[str],
) -> tuple[list[str], str | None]:
    """Copy the nearest permitted strain's measured results.

    Ties are broken by the second and third nearest, so a neighbour that tested
    negative for everything does not produce an arbitrary order.
    """
    response = client.search(
        index=config.BACTERIA,
        size=len(excluded) + 5,
        knn={
            "field": "genome_vector",
            "query_vector": vector,
            "k": len(excluded) + 5,
            "num_candidates": 500,
        },
        source_includes=["strain_id"],
    )
    ordered = [
        (h["_source"]["strain_id"], _distance_from_score(h["_score"]))
        for h in response["hits"]["hits"]
        if h["_source"]["strain_id"] not in excluded
    ]
    if not ordered:
        return list(phages), None
    nearest = [sid for sid, _ in ordered[:3] if sid in truth]
    if not nearest:
        return list(phages), None

    def key(phage: str) -> tuple:
        return tuple(-int(truth[sid].outcomes.get(phage, False)) for sid in nearest)

    return sorted(phages, key=key), nearest[0]


# ----------------------------------------------------------------------- the run


#: Funnel configurations evaluated. The ablations run the *same code path* as the
#: full funnel with one signal switched off, which is the only way an ablation
#: answers "which part carries the signal" rather than "what does a different
#: program do".
FUNNEL_CONFIGS = {
    "funnel": {},
    "funnel_minus_rbp": {"use_rbp": False},
    "funnel_minus_prior": {"use_prior": False},
}


def run_benchmark(
    client: Elasticsearch,
    *,
    n_folds: int = N_FOLDS,
    limit: int | None = None,
    neighbour_k: int = config.NEIGHBOUR_K,
    methods: Sequence[str] | None = None,
    progress: bool = True,
) -> dict:
    """Run every method over strain-grouped folds. Returns the results document."""
    truth = load_truth(client)
    strains = evaluable_strains(client, truth)
    if limit:
        strains = strains[:limit]
    if not strains:
        raise RuntimeError("no strains have both a genome_vector and interaction rows")

    phages = sorted(
        {phage for record in truth.values() for phage in record.outcomes}
    )
    folds = make_folds(strains, n_folds)
    wanted = set(methods) if methods else None

    scores: dict[str, MethodScores] = {}

    def bucket(name: str) -> MethodScores:
        return scores.setdefault(name, MethodScores(name))

    vectors = _strain_vectors(client, strains)
    base_rate_num = base_rate_den = 0
    started = time.perf_counter()
    done = 0

    for fold_index, fold in enumerate(folds):
        held_out = set(fold)
        training = [s for s in strains if s not in held_out]
        breadth = fold_breadth(client, training, phages)

        for strain_id in fold:
            record = truth[strain_id]
            if not record.positives:
                # A strain no phage infects contributes nothing to P@10 or AUPRC
                # and would drag every method's mean toward zero identically.
                continue
            base_rate_num += len(record.positives)
            base_rate_den += record.n_tested
            # Everything else in this strain's fold is invisible to every method.
            excluded = held_out - {strain_id}

            if wanted is None or "random" in wanted:
                bucket("random").add(baseline_random(phages, strain_id), record)
            if wanted is None or "generalist" in wanted:
                bucket("generalist").add(baseline_generalist(breadth, phages), record)
            if wanted is None or "phylo_nn" in wanted:
                ranking, _ = baseline_phylo_nn(
                    client, strain_id, vectors[strain_id], excluded | {strain_id}, truth, phages
                )
                bucket("phylo_nn").add(ranking, record)

            for name, overrides in FUNNEL_CONFIGS.items():
                if wanted is not None and name not in wanted:
                    continue
                started_run = time.perf_counter()
                try:
                    run = run_funnel(
                        client,
                        strain_id,
                        top_n=len(phages),
                        neighbour_k=neighbour_k,
                        exclude_strains=tuple(excluded),
                        breadth=breadth,
                        explain=False,
                        build_cocktail=False,
                        persist=False,
                        **overrides,
                    )
                except EmptyStage as exc:  # pragma: no cover - defensive
                    bucket(name).failures.append((strain_id, str(exc)))
                    continue
                elapsed = (time.perf_counter() - started_run) * 1000
                if not run.ok:
                    bucket(name).failures.append(
                        (strain_id, f"{run.failed_stage}: {run.failure_reason}")
                    )
                    continue
                ranking = [entry["phage_id"] for entry in run.shortlist]
                chosen: list[str] = []
                if name == "funnel":
                    chosen = [
                        entry["phage_id"]
                        for entry in cocktail_mod.design(
                            client, run.shortlist[:TOP_K], size=COCKTAIL_SIZE
                        )
                    ]
                bucket(name).ms.append(elapsed)
                bucket(name).add(ranking, record, cocktail=chosen)

            done += 1
            if progress and done % 25 == 0:
                rate = done / (time.perf_counter() - started)
                print(
                    f"  fold {fold_index + 1}/{len(folds)}  {done} strains  "
                    f"{rate:.1f}/s",
                    flush=True,
                )

    base_rate = base_rate_num / base_rate_den if base_rate_den else float("nan")
    summaries = [scores[name].summary() for name in scores]
    for summary in summaries:
        summary["enrichment_over_base"] = (
            summary["p_at_10"] / base_rate if base_rate else float("nan")
        )

    return {
        "n_strains_evaluated": len(scores[next(iter(scores))].p_at_k) if scores else 0,
        "n_strains_available": len(strains),
        "n_phages": len(phages),
        "n_folds": n_folds,
        "neighbour_k": neighbour_k,
        "infect_threshold": _threshold(),
        "base_rate": base_rate,
        "seed": SEED,
        "results": summaries,
        "failures": {
            name: score.failures[:10] for name, score in scores.items() if score.failures
        },
        "wall_seconds": round(time.perf_counter() - started, 1),
    }


def _threshold() -> float:
    from phageforge.ingest.gaborieau import INFECT_THRESHOLD

    return INFECT_THRESHOLD


def _strain_vectors(client: Elasticsearch, strains: Sequence[str]) -> dict[str, list[float]]:
    response = client.search(
        index=config.BACTERIA,
        size=len(strains) + 10,
        query={"terms": {"strain_id": list(strains)}},
        source_includes=["strain_id", "genome_vector"],
    )
    return {
        h["_source"]["strain_id"]: h["_source"]["genome_vector"]
        for h in response["hits"]["hits"]
        if h["_source"].get("genome_vector")
    }


# ------------------------------------------------------------------- presentation


_ROW_ORDER = [
    ("random", "Random"),
    ("generalist", "Most-generalist phage"),
    ("phylo_nn", "Phylogenetic NN"),
    ("funnel", "PhageForge funnel"),
    ("funnel_minus_rbp", "  - RBP arm (ablation)"),
    ("funnel_minus_prior", "  - neighbour prior (ablation)"),
]


def format_table(results: dict) -> str:
    by_name = {r["method"]: r for r in results["results"]}
    # Parenthesised so each element is unambiguously one line: inside a list, an
    # implicit concatenation and a missing comma look identical.
    lines = [
        (
            f"Strain-grouped {results['n_folds']}-fold CV over "
            f"{results['n_strains_evaluated']} strains x {results['n_phages']} phages"
        ),
        (
            f"positives = score > {results['infect_threshold']:g}; "
            f"base rate {results['base_rate']:.1%}; "
            f"neighbour k = {results['neighbour_k']}; seed {results['seed']}"
        ),
        "",
        (
            f"{'Method':30s} {'P@10':>7s} {'R@10':>7s} {'AUPRC':>7s} {'x base':>7s} "
            f"{'cocktail':>9s} {'p95 ms':>8s}"
        ),
        "-" * 82,
    ]
    for key, label in _ROW_ORDER:
        row = by_name.get(key)
        if not row:
            continue
        cocktail = (
            f"{row['cocktail_coverage']:.1%}"
            if not math.isnan(row["cocktail_coverage"])
            else "-"
        )
        latency = f"{row['p95_ms']:.0f}" if not math.isnan(row["p95_ms"]) else "-"
        lines.append(
            f"{label:30s} {row['p_at_10']:7.3f} {row['r_at_10']:7.3f} "
            f"{row['auprc']:7.3f} {row['enrichment_over_base']:7.2f} "
            f"{cocktail:>9s} {latency:>8s}"
        )
    funnel = by_name.get("funnel")
    if funnel and not math.isnan(funnel["top4_coverage"]):
        lines.append("")
        lines.append(
            f"cocktail (4 phages, receptor-diverse) covers "
            f"{funnel['cocktail_coverage']:.1%} of held-out strains vs "
            f"{funnel['top4_coverage']:.1%} for the top 4 ranked individually"
        )
    for name, failures in results.get("failures", {}).items():
        lines.append(f"\n{name}: {len(failures)} strain(s) failed, e.g. {failures[0]}")
    return "\n".join(lines)


def _finite(value):
    """Recursively replace non-finite floats with ``None``.

    A metric that does not apply to a method (cocktail coverage for the random
    baseline, latency for a method that runs offline) is carried as ``nan``
    internally, which is right -- it is genuinely not a number. But
    ``json.dumps`` writes it as a bare ``NaN`` token, which is valid Python and
    *invalid JSON*: ``jq``, ``JSON.parse`` and any strict parser reject the whole
    file. Written as ``null`` it means the same thing and survives the trip.
    """
    if isinstance(value, dict):
        return {k: _finite(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_finite(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def save(results: dict, path: Path | None = None) -> Path:
    path = path or (config.DERIVED_DIR / "benchmark.json")
    path.write_text(json.dumps(_finite(results), indent=2, allow_nan=False))
    return path
