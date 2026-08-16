"""Response schemas for the PhageForge API.

These mirror what the funnel already produces rather than reshaping it. Two
choices are load-bearing and deliberate:

* Every shortlist entry carries ``reasons`` and an ``evidence_tier``. A
  recommendation without a stated basis is not a deliverable product
  (TECHNICAL_DESIGN.md G3), so the field is required, not optional.
* Nothing here is called a probability. Stage D emits an uncalibrated linear
  score and the cocktail uses a ``probability_proxy``; the names say so, because
  a number labelled "probability" will be read as one.
"""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field


class Taxonomy(BaseModel):
    realm: str | None = None
    family: str | None = None
    subfamily: str | None = None
    genus: str | None = None
    species: str | None = None


class Signal(BaseModel):
    """One significant_terms hit from Stage C."""

    field: str
    value: str
    fg_count: int = Field(description="hosts of this phage carrying the feature")
    fg_total: int = Field(description="hosts of this phage that were tested")
    bg_count: int = Field(description="all strains carrying the feature")
    bg_total: int
    score: float = Field(description="Elasticsearch significance score, not a p-value")
    matches_target: bool = Field(description="the queried strain carries this feature")


class ShortlistEntry(BaseModel):
    rank: int
    phage_id: str
    score: float = Field(description="uncalibrated Stage D linear score; ranking only")
    evidence_tier: int = Field(
        description="1 = measured strain-level interactions among neighbours; "
        "4 = vector similarity alone"
    )
    prior: float = Field(description="similarity-weighted neighbour infection rate, 0-1")
    n_neighbours_tested: int
    n_neighbours_infected: int
    taxonomy: Taxonomy = Taxonomy()
    morphotype: str | None = None
    lifestyle: str | None = None
    genome_length: int | None = None
    host_range_breadth: float | None = None
    rbp_confidence: str | None = Field(
        default=None, description="high | heuristic — how the RBP locus tag resolved"
    )
    reasons: list[str] = Field(description="human-readable basis; never empty")
    signals: list[Signal] = []
    features: dict[str, float] = Field(
        default_factory=dict, description="the Stage D feature vector, for audit"
    )


class CocktailEntry(BaseModel):
    phage_id: str
    probability_proxy: float = Field(
        description="neighbour prior, or normalised score when untested. Not calibrated."
    )
    receptor_overlap: float = Field(description="summed RBP cosine against the chosen set")
    marginal_gain: float
    rbp_confidence: str | None = None
    overlap_basis: str = Field(description="esm2_rbp_cosine | none")
    rationale: str


class StageTiming(BaseModel):
    ms: float
    n_results: int
    metadata: dict[str, Any] = {}


class StrainSummary(BaseModel):
    strain_id: str
    species: str | None = None
    st: str | None = None
    phylogroup: str | None = None
    o_antigen: str | None = None
    h_type: str | None = None
    lps_type: str | None = None
    pathotype: str | None = None
    n_infections: int | None = None
    defence_systems: list[str] = []
    has_genome_vector: bool = True


class PredictionResponse(BaseModel):
    run_id: str
    strain_id: str
    strain: StrainSummary
    model_version: str
    shortlist: list[ShortlistEntry]
    cocktail: list[CocktailEntry] = []
    stages: dict[str, StageTiming] = {}
    total_ms: float
    persisted: bool = Field(description="whether the run was written to pf-predictions")


class FunnelStopped(BaseModel):
    """A funnel that produced nothing, with the stage and cause named.

    An empty shortlist is never returned as success. TECHNICAL_DESIGN.md's rule
    is that a stage yielding nothing raises a named condition rather than a blank
    list, and this is that condition crossing the HTTP boundary.
    """

    detail: str
    strain_id: str
    failed_stage: str
    failure_reason: str
    run_id: str


class PhageSummary(BaseModel):
    phage_id: str
    name: str | None = None
    accession: str | None = None
    taxonomy: Taxonomy = Taxonomy()
    morphotype: str | None = None
    lifestyle: str | None = None
    genome_length: int | None = None
    host_range_breadth: float | None = None
    n_rbps: int = 0
    rbp_match_confidence: str | None = None
    has_rbp_vector: bool = False


class IndexStatus(BaseModel):
    docs: int
    exists: bool


class HealthResponse(BaseModel):
    status: str = Field(description="ok | degraded")
    cluster: str
    cluster_status: str
    model_version: str
    indices: dict[str, IndexStatus]
    ready: bool = Field(description="the funnel has everything it needs to answer")
    notes: list[str] = []
