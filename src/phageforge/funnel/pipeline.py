"""Funnel orchestration: run the stages, persist the run, return the shortlist.

Thin by design. This module decides *what runs in what order* and records what
happened; it does no scoring of its own. The one exception is the cocktail
set-cover, which is a combinatorial optimisation over the finished shortlist and
is the single carve-out TECHNICAL_DESIGN.md R3 allows.
"""

from __future__ import annotations

import uuid
from dataclasses import dataclass, field
from datetime import UTC, datetime

from elasticsearch import Elasticsearch

from phageforge import config
from phageforge.es import bulk
from phageforge.funnel import cocktail as cocktail_mod
from phageforge.funnel.stages import (
    MODEL_VERSION,
    EmptyStage,
    Prior,
    StageResult,
    find_neighbours,
    get_strain,
    stage_a_candidates,
    stage_b_prior,
    stage_c_explain,
    stage_d_rank,
)


@dataclass
class FunnelRun:
    run_id: str
    strain_id: str
    strain: dict
    shortlist: list[dict]
    stages: dict[str, StageResult] = field(default_factory=dict)
    cocktail: list[dict] = field(default_factory=list)
    total_ms: float = 0.0
    failed_stage: str | None = None
    failure_reason: str | None = None

    @property
    def ok(self) -> bool:
        return self.failed_stage is None

    def stage_ms(self) -> dict[str, float]:
        return {name.lower(): round(r.ms, 2) for name, r in self.stages.items()}

    def to_document(self) -> dict:
        """The ``pf-predictions`` audit record for this run."""
        return {
            "run_id": self.run_id,
            "strain_id": self.strain_id,
            "created_at": datetime.now(UTC).isoformat(),
            "model_version": MODEL_VERSION,
            "stages": {
                name: {"ms": round(r.ms, 2), **r.metadata} for name, r in self.stages.items()
            },
            "shortlist": self.shortlist,
            "shortlist_ids": [p["phage_id"] for p in self.shortlist],
            "n_candidates": len(self.stages["A"].results) if "A" in self.stages else 0,
            "n_shortlist": len(self.shortlist),
            "stage_ms": self.stage_ms(),
            "total_ms": round(self.total_ms, 2),
            "failed_stage": self.failed_stage,
            "failure_reason": self.failure_reason,
        }


def run_funnel(
    client: Elasticsearch,
    strain_id: str,
    *,
    top_n: int = config.TOP_N,
    neighbour_k: int = config.NEIGHBOUR_K,
    candidate_pool: int = config.CANDIDATE_POOL,
    exclude_strains: tuple[str, ...] = (),
    use_rbp: bool = True,
    use_prior: bool = True,
    weights: dict[str, float] | None = None,
    breadth: dict[str, float] | None = None,
    explain: bool = True,
    build_cocktail: bool = True,
    persist: bool = True,
) -> FunnelRun:
    """Run N -> B -> A -> C -> D for one strain.

    ``use_rbp`` and ``use_prior`` are the two ablation switches from
    TECHNICAL_DESIGN.md 7.6, wired here rather than bolted onto the benchmark so
    that the ablated funnel is the *same code path* as the real one.

    ``exclude_strains`` is how the benchmark prevents leakage: the held-out fold
    is excluded from the neighbour set, so a strain is never scored using its own
    fold's measurements.

    A stage that produces nothing raises :class:`EmptyStage`; this function
    catches it, records which stage failed and why, and returns a run with an
    empty shortlist and ``ok == False``. The caller must never see a blank
    shortlist without a stated cause.
    """
    run = FunnelRun(run_id=uuid.uuid4().hex[:16], strain_id=strain_id, strain={}, shortlist=[])
    started = datetime.now(UTC)

    try:
        strain = get_strain(client, strain_id)
        run.strain = {k: v for k, v in strain.items() if k != "genome_vector"}

        neighbours = find_neighbours(client, strain, k=neighbour_k, exclude=exclude_strains)
        run.stages["N"] = neighbours

        prior = stage_b_prior(client, neighbours.results)
        run.stages["B"] = prior
        priors: list[Prior] = prior.results
        if not use_prior:
            # Ablation: keep the phage universe Stage B discovered, drop its
            # signal. Zeroing rather than skipping keeps Stage A seedable, so the
            # ablation isolates the prior instead of collapsing the whole funnel.
            priors = [Prior(p.phage_id, 0.0, p.hit_rate, p.n_tested, p.n_hits, p.weight_sum)
                      for p in priors]

        candidates = stage_a_candidates(
            client,
            prior.results,  # seeding always uses the true prior; the ablation is on scoring
            pool=candidate_pool,
            use_rbp=use_rbp,
        )
        run.stages["A"] = candidates

        ranked = stage_d_rank(
            client, candidates.results, priors, top_n=top_n, weights=weights, breadth=breadth
        )
        run.stages["D"] = ranked
        run.shortlist = ranked.results

        if explain:
            explanations = stage_c_explain(
                client, [p["phage_id"] for p in run.shortlist], strain
            )
            run.stages["C"] = explanations
            for entry in run.shortlist:
                detail = explanations.results.get(entry["phage_id"], {})
                entry["reasons"] = detail.get("reasons", [])
                entry["signals"] = detail.get("signals", [])[:5]
                if not entry["reasons"]:
                    # G3 requires a reason on every result. When no feature of the
                    # target is enriched among this phage's hosts, say that
                    # plainly rather than emitting an empty list.
                    entry["reasons"] = [_fallback_reason(entry)]

        if build_cocktail and run.shortlist:
            run.cocktail = cocktail_mod.design(client, run.shortlist)

    except EmptyStage as exc:
        run.failed_stage = exc.stage
        run.failure_reason = exc.reason

    run.total_ms = (datetime.now(UTC) - started).total_seconds() * 1000

    if persist:
        bulk.bulk_index(
            client,
            config.PREDICTIONS,
            [(bulk.doc_id("run", run.run_id), run.to_document())],
        )
    return run


def _fallback_reason(entry: dict) -> str:
    if entry["n_neighbours_tested"]:
        return (
            f"no enriched feature of this strain among its known hosts; ranked on "
            f"neighbour evidence ({entry['n_neighbours_infected']}/"
            f"{entry['n_neighbours_tested']} genomically similar strains infected)"
        )
    return (
        "no measured evidence in this strain's neighbourhood; ranked on genome and "
        "RBP similarity alone (evidence tier 4)"
    )


def format_run(run: FunnelRun) -> str:
    """Human-readable funnel report for the CLI."""
    lines: list[str] = []
    strain = run.strain
    lines.append(f"strain {run.strain_id}  run {run.run_id}")
    lines.append(
        "  "
        + "  ".join(
            f"{k}={strain.get(k)}"
            for k in ("st", "phylogroup", "o_antigen", "lps_type")
            if strain.get(k)
        )
    )

    if not run.ok:
        lines.append(f"\n  FUNNEL STOPPED at stage {run.failed_stage}: {run.failure_reason}")
        return "\n".join(lines)

    counts = {
        "N": len(run.stages["N"].results),
        "B": len(run.stages["B"].results),
        "A": len(run.stages["A"].results),
        "D": len(run.stages["D"].results),
    }
    lines.append("")
    lines.append(
        f"  N neighbours  {counts['N']:>5}  ({run.stages['N'].ms:.0f} ms)  "
        f"distance {run.stages['N'].metadata['distance_min']:.4f}"
        f"-{run.stages['N'].metadata['distance_max']:.4f}, "
        f"{run.stages['N'].metadata['same_st']} share the target ST"
    )
    lines.append(
        f"  B prior       {counts['B']:>5}  ({run.stages['B'].ms:.0f} ms)  "
        f"{run.stages['B'].metadata['n_with_positive_prior']} with a positive prior"
    )
    lines.append(
        f"  A candidates  {counts['A']:>5}  ({run.stages['A'].ms:.0f} ms)  "
        f"arms: {'+'.join(run.stages['A'].metadata['arms'])}, "
        f"corpus {run.stages['A'].metadata['corpus_size']}"
    )
    if "C" in run.stages:
        lines.append(
            f"  C explanation {run.stages['C'].metadata['n_phages']:>5}  "
            f"({run.stages['C'].ms:.0f} ms)  "
            f"{run.stages['C'].metadata['n_with_matching_reason']} with a feature "
            "matching the target"
        )
    lines.append(
        f"  D ranking     {counts['D']:>5}  ({run.stages['D'].ms:.0f} ms)  "
        f"{run.stages['D'].metadata['model_version']}"
    )
    lines.append(f"  total {run.total_ms:.0f} ms")

    lines.append("\n  shortlist:")
    for i, entry in enumerate(run.shortlist, start=1):
        genus = (entry.get("taxonomy") or {}).get("genus") or "?"
        lines.append(
            f"   {i:2d}. {entry['phage_id']:14s} score {entry['score']:6.3f}  "
            f"tier {entry['evidence_tier']}  prior {entry['prior']:.2f}  "
            f"({entry['n_neighbours_infected']}/{entry['n_neighbours_tested']} neighbours)  "
            f"{genus}"
        )
        for reason in entry.get("reasons", [])[:2]:
            lines.append(f"       - {reason}")

    if run.cocktail:
        lines.append("\n  cocktail:")
        for entry in run.cocktail:
            lines.append(
                f"   + {entry['phage_id']:14s} marginal gain {entry['marginal_gain']:.3f}  "
                f"{entry['rationale']}"
            )
    return "\n".join(lines)
