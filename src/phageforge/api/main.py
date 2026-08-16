"""FastAPI service over the funnel (Phase 5).

Thin on purpose. The funnel already decides what runs, records what happened and
refuses to return an unexplained shortlist; this layer adds HTTP, validation and
lookup endpoints, and does no biology of its own. If a behaviour matters, it
belongs in ``funnel/`` where the CLI and the benchmark get it too.

Two rules this layer does enforce, because they are HTTP-shaped:

* **An empty shortlist is never a 200.** ``EmptyStage`` names the stage and the
  cause; that becomes a 422 carrying both, so a caller can never mistake "the
  funnel stopped" for "no phage is suitable".
* **An unknown strain is a 404, not a funnel failure.** The distinction matters
  to a client deciding whether to retry.
"""

from __future__ import annotations

import json
import math
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any

from elasticsearch import Elasticsearch
from fastapi import Depends, FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse

STATIC_DIR = Path(__file__).parent / "static"

from phageforge import config
from phageforge.api import models
from phageforge.es import bulk
from phageforge.funnel import pipeline
from phageforge.funnel.stages import MODEL_VERSION

_state: dict[str, Any] = {}


@asynccontextmanager
async def lifespan(app: FastAPI):
    """One Elasticsearch client for the process lifetime.

    Health is checked but not required at startup: a service that refuses to
    boot when the cluster is briefly unavailable cannot report *why* it is
    unavailable. ``/health`` tells the truth instead.
    """
    client = config.get_client()
    _state["client"] = client
    try:
        _state["health"] = config.require_healthy(client)
    except Exception as exc:  # noqa: BLE001 - reported through /health, not raised
        _state["health"] = None
        _state["health_error"] = str(exc)
    yield
    client.close()


app = FastAPI(
    title="PhageForge",
    version="0.1.0",
    summary="Hybrid retrieval funnel for phage-host prediction, built on Elasticsearch",
    description=(
        "Given an *E. coli* strain, return a ranked, explained shortlist of phages "
        "likely to infect it, plus a receptor-diverse cocktail.\n\n"
        "Scores are uncalibrated ranking scores, not probabilities. Every result "
        "carries an evidence tier and a stated reason."
    ),
    lifespan=lifespan,
)


def get_es() -> Elasticsearch:
    client = _state.get("client")
    if client is None:  # pragma: no cover - only before lifespan runs
        raise HTTPException(status_code=503, detail="Elasticsearch client not initialised")
    return client


# ----------------------------------------------------------------- helpers


def _strain_doc(client: Elasticsearch, strain_id: str) -> dict:
    """Fetch one strain by id, or raise 404.

    Done here rather than leaving it to the funnel so that "no such strain" is a
    404 with a suggestion, not a 422 about stage N.
    """
    hits = client.search(
        index=config.BACTERIA,
        size=1,
        query={"term": {"strain_id": strain_id}},
        source_excludes=["genome_vector"],
    )["hits"]["hits"]
    if not hits:
        near = _suggest_strains(client, strain_id)
        detail = f"strain {strain_id!r} is not in {config.BACTERIA}"
        if near:
            detail += f". Did you mean: {', '.join(near)}?"
        raise HTTPException(status_code=404, detail=detail)
    return hits[0]["_source"]


def _suggest_strains(client: Elasticsearch, strain_id: str, limit: int = 5) -> list[str]:
    """Cheap did-you-mean over strain ids. Best effort; never raises."""
    try:
        hits = client.search(
            index=config.BACTERIA,
            size=limit,
            query={"wildcard": {"strain_id": {"value": f"*{strain_id}*", "case_insensitive": True}}},
            source_includes=["strain_id"],
        )["hits"]["hits"]
    except Exception:  # noqa: BLE001 - a suggestion failing must not mask the 404
        return []
    return [h["_source"]["strain_id"] for h in hits]


def _strain_summary(source: dict, *, has_vector: bool | None = None) -> models.StrainSummary:
    return models.StrainSummary(
        strain_id=source["strain_id"],
        species=source.get("species"),
        st=source.get("st"),
        phylogroup=source.get("phylogroup"),
        o_antigen=source.get("o_antigen"),
        h_type=source.get("h_type"),
        lps_type=source.get("lps_type"),
        pathotype=source.get("pathotype"),
        n_infections=source.get("n_infections"),
        defence_systems=source.get("defence_systems") or [],
        has_genome_vector=(
            has_vector if has_vector is not None else bool(source.get("genome_vector"))
        ),
    )


def _shortlist(run: pipeline.FunnelRun) -> list[models.ShortlistEntry]:
    return [
        models.ShortlistEntry(
            rank=i,
            phage_id=e["phage_id"],
            score=e["score"],
            evidence_tier=e["evidence_tier"],
            prior=e["prior"],
            n_neighbours_tested=e["n_neighbours_tested"],
            n_neighbours_infected=e["n_neighbours_infected"],
            taxonomy=models.Taxonomy(**(e.get("taxonomy") or {})),
            morphotype=e.get("morphotype"),
            lifestyle=e.get("lifestyle"),
            genome_length=e.get("genome_length"),
            host_range_breadth=e.get("host_range_breadth"),
            rbp_confidence=e.get("rbp_confidence"),
            reasons=e.get("reasons") or [],
            signals=[models.Signal(**_signal(s)) for s in e.get("signals", [])],
            features=e.get("features") or {},
        )
        for i, e in enumerate(run.shortlist, start=1)
    ]


def _signal(signal: dict) -> dict:
    keys = ("field", "value", "fg_count", "fg_total", "bg_count", "bg_total",
            "score", "matches_target")
    return {k: signal[k] for k in keys}


def _stage_timings(run: pipeline.FunnelRun) -> dict[str, models.StageTiming]:
    return {
        name: models.StageTiming(
            ms=round(result.ms, 2),
            n_results=len(result.results),
            metadata=_jsonable(result.metadata),
        )
        for name, result in run.stages.items()
    }


def _jsonable(value: Any) -> Any:
    """Strip non-finite floats. ``NaN`` is valid Python and invalid JSON."""
    if isinstance(value, dict):
        return {k: _jsonable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


# ----------------------------------------------------------------- endpoints


@app.get("/health", response_model=models.HealthResponse, tags=["service"])
def health(client: Elasticsearch = Depends(get_es)) -> models.HealthResponse:
    """Cluster reachability, per-index counts, and whether the funnel can answer.

    ``ready`` is the useful field: the cluster can be green while the funnel is
    unable to serve a single request because features were never computed.
    """
    notes: list[str] = []
    try:
        cluster = dict(client.cluster.health())
    except Exception as exc:  # noqa: BLE001
        return models.HealthResponse(
            status="degraded",
            cluster="unreachable",
            cluster_status="unknown",
            model_version=MODEL_VERSION,
            indices={},
            ready=False,
            notes=[f"Elasticsearch unreachable: {exc}"],
        )

    indices: dict[str, models.IndexStatus] = {}
    for index in config.LOCAL_INDICES:
        exists = client.indices.exists(index=index)
        indices[index] = models.IndexStatus(
            docs=bulk.count(client, index) if exists else 0, exists=bool(exists)
        )

    vectors = 0
    if indices[config.BACTERIA].exists:
        vectors = bulk.count(client, config.BACTERIA, {"exists": {"field": "genome_vector"}})
        missing = indices[config.BACTERIA].docs - vectors
        if missing:
            notes.append(
                f"{missing} strain{'s' if missing > 1 else ''} "
                f"{'have' if missing > 1 else 'has'} no genome_vector and cannot be queried"
            )
    if not indices[config.INTERACTIONS].docs:
        notes.append("pf-interactions is empty; the Stage B prior has nothing to transfer")
    if not indices[config.PROTEINS].docs:
        notes.append("pf-proteins is empty; Stage A runs without its RBP arm")

    ready = bool(vectors and indices[config.INTERACTIONS].docs and indices[config.PHAGES].docs)
    return models.HealthResponse(
        status="ok" if ready and cluster["status"] != "red" else "degraded",
        cluster=cluster["cluster_name"],
        cluster_status=cluster["status"],
        model_version=MODEL_VERSION,
        indices=indices,
        ready=ready,
        notes=notes,
    )


@app.get("/strains", response_model=list[models.StrainSummary], tags=["catalogue"])
def list_strains(
    q: str | None = Query(default=None, description="substring of the strain id"),
    phylogroup: str | None = None,
    st: str | None = None,
    limit: int = Query(default=25, ge=1, le=500),
    client: Elasticsearch = Depends(get_es),
) -> list[models.StrainSummary]:
    """Browse or type-ahead the strain catalogue."""
    filters: list[dict] = []
    if q:
        filters.append({"wildcard": {"strain_id": {"value": f"*{q}*", "case_insensitive": True}}})
    if phylogroup:
        filters.append({"term": {"phylogroup": phylogroup}})
    if st:
        filters.append({"term": {"st": st}})

    hits = client.search(
        index=config.BACTERIA,
        size=limit,
        query={"bool": {"filter": filters}} if filters else {"match_all": {}},
        sort=[{"strain_id": "asc"}],
        source_excludes=["genome_vector"],
    )["hits"]["hits"]
    return [_strain_summary(h["_source"], has_vector=True) for h in hits]


@app.get("/strains/{strain_id}", response_model=models.StrainSummary, tags=["catalogue"])
def get_strain(strain_id: str, client: Elasticsearch = Depends(get_es)) -> models.StrainSummary:
    source = _strain_doc(client, strain_id)
    # source_excludes drops genome_vector, so ask the index directly.
    has_vector = bool(
        client.search(
            index=config.BACTERIA,
            size=0,
            query={
                "bool": {
                    "filter": [
                        {"term": {"strain_id": strain_id}},
                        {"exists": {"field": "genome_vector"}},
                    ]
                }
            },
        )["hits"]["total"]["value"]
    )
    return _strain_summary(source, has_vector=has_vector)


@app.get("/phages", response_model=list[models.PhageSummary], tags=["catalogue"])
def list_phages(
    q: str | None = Query(default=None, description="substring of the phage id"),
    genus: str | None = None,
    limit: int = Query(default=25, ge=1, le=500),
    client: Elasticsearch = Depends(get_es),
) -> list[models.PhageSummary]:
    filters: list[dict] = []
    if q:
        filters.append({"wildcard": {"phage_id": {"value": f"*{q}*", "case_insensitive": True}}})
    if genus:
        filters.append({"term": {"taxonomy.genus": genus}})

    hits = client.search(
        index=config.PHAGES,
        size=limit,
        query={"bool": {"filter": filters}} if filters else {"match_all": {}},
        sort=[{"phage_id": "asc"}],
        source_excludes=["genome_vector", "rbp_match_vector"],
    )["hits"]["hits"]
    return [_phage_summary(h["_source"]) for h in hits]


@app.get("/phages/{phage_id}", response_model=models.PhageSummary, tags=["catalogue"])
def get_phage(phage_id: str, client: Elasticsearch = Depends(get_es)) -> models.PhageSummary:
    hits = client.search(
        index=config.PHAGES,
        size=1,
        query={"term": {"phage_id": phage_id}},
        source_excludes=["genome_vector", "rbp_match_vector"],
    )["hits"]["hits"]
    if not hits:
        raise HTTPException(status_code=404, detail=f"phage {phage_id!r} is not in {config.PHAGES}")
    return _phage_summary(hits[0]["_source"])


def _phage_summary(source: dict) -> models.PhageSummary:
    return models.PhageSummary(
        phage_id=source["phage_id"],
        name=source.get("name"),
        accession=source.get("accession"),
        taxonomy=models.Taxonomy(**(source.get("taxonomy") or {})),
        morphotype=source.get("morphotype"),
        lifestyle=source.get("lifestyle"),
        genome_length=source.get("genome_length"),
        host_range_breadth=source.get("host_range_breadth"),
        n_rbps=len(source.get("rbp_ids") or []),
        rbp_match_confidence=source.get("rbp_match_confidence"),
        # rbp_match_vector is excluded from _source above, so infer presence from
        # the id recorded alongside it (ES 9 omits dense_vector by default anyway).
        has_rbp_vector=bool(source.get("rbp_match_id")),
    )


@app.get(
    "/predict/{strain_id}",
    response_model=models.PredictionResponse,
    responses={404: {"description": "unknown strain"},
               422: {"model": models.FunnelStopped, "description": "the funnel stopped"}},
    tags=["predict"],
)
def predict(
    strain_id: str,
    top_n: int = Query(default=config.TOP_N, ge=1, le=96),
    neighbours: int = Query(default=config.NEIGHBOUR_K, ge=1, le=402,
                            description="k for the Stage N kNN over strain genome vectors"),
    explain: bool = Query(default=True, description="run Stage C for reason strings"),
    cocktail: bool = Query(default=True, description="run the greedy set-cover"),
    use_rbp: bool = Query(default=True, description="ablation: Stage A's RBP arm"),
    use_prior: bool = Query(default=True, description="ablation: the Stage B prior"),
    persist: bool = Query(default=True, description="write the run to pf-predictions"),
    client: Elasticsearch = Depends(get_es),
):
    """Rank phages for one strain.

    The ablation switches are exposed because they are the same code path the
    benchmark measures — a reviewer can reproduce an ablation row over HTTP
    rather than taking the benchmark's word for it.
    """
    _strain_doc(client, strain_id)  # 404 before the funnel, so the error is honest

    run = pipeline.run_funnel(
        client,
        strain_id,
        top_n=top_n,
        neighbour_k=neighbours,
        use_rbp=use_rbp,
        use_prior=use_prior,
        explain=explain,
        build_cocktail=cocktail,
        persist=persist,
    )

    if not run.ok:
        return JSONResponse(
            status_code=422,
            content=models.FunnelStopped(
                detail=f"the funnel stopped at stage {run.failed_stage}",
                strain_id=strain_id,
                failed_stage=run.failed_stage or "?",
                failure_reason=run.failure_reason or "unknown",
                run_id=run.run_id,
            ).model_dump(),
        )

    return models.PredictionResponse(
        run_id=run.run_id,
        strain_id=run.strain_id,
        strain=_strain_summary(run.strain, has_vector=True),
        model_version=MODEL_VERSION,
        shortlist=_shortlist(run),
        cocktail=[models.CocktailEntry(**c) for c in run.cocktail],
        stages=_stage_timings(run),
        total_ms=round(run.total_ms, 2),
        persisted=persist,
    )


@app.get("/runs", tags=["audit"])
def list_runs(
    strain_id: str | None = None,
    limit: int = Query(default=20, ge=1, le=200),
    client: Elasticsearch = Depends(get_es),
) -> list[dict]:
    """Recent funnel runs, newest first. The audit trail, summary form."""
    query = {"term": {"strain_id": strain_id}} if strain_id else {"match_all": {}}
    hits = client.search(
        index=config.PREDICTIONS,
        size=limit,
        query=query,
        sort=[{"created_at": "desc"}],
        source_includes=["run_id", "strain_id", "created_at", "model_version",
                         "n_shortlist", "shortlist_ids", "total_ms", "failed_stage"],
    )["hits"]["hits"]
    return [_jsonable(h["_source"]) for h in hits]


@app.get("/runs/{run_id}", tags=["audit"])
def get_run(run_id: str, client: Elasticsearch = Depends(get_es)) -> dict:
    """One stored run in full, including per-stage metadata."""
    hits = client.search(
        index=config.PREDICTIONS, size=1, query={"term": {"run_id": run_id}}
    )["hits"]["hits"]
    if not hits:
        raise HTTPException(status_code=404, detail=f"run {run_id!r} not found")
    return _jsonable(hits[0]["_source"])


@app.get("/benchmark", tags=["audit"])
def benchmark(
    threshold: int = Query(
        default=0,
        ge=0,
        le=3,
        description="which infection-threshold run to serve. The code binarises at "
        "`score > threshold`, so 0 = any lysis, 1 = the '>=2' cut, 2 = the '>=3' cut.",
    ),
) -> dict:
    """A measured benchmark run, as written by ``make bench``.

    Served from disk rather than recomputed: the benchmark takes ~70 s and its
    whole point is that the number you are shown is the number that was measured.

    Each threshold is a *separate artefact* because changing it requires a full
    re-derivation (``ingest -> features -> proteins -> bench``) — ``infects`` is
    computed at ingest, so the runs cannot be produced by re-reading one file.
    """
    name = "benchmark.json" if threshold == 0 else f"benchmark-gt{threshold}.json"
    path = config.DERIVED_DIR / name
    if not path.exists():
        raise HTTPException(
            status_code=404,
            detail=(
                f"no {name} yet — run `PF_INFECT_THRESHOLD={threshold} make ingest features "
                f"proteins bench` (see current_status.md §3.1)"
                if threshold
                else "no benchmark.json yet — run `make bench`"
            ),
        )
    # parse_constant catches the bare NaN that older runs wrote before the
    # harness switched to emitting null.
    results = json.loads(path.read_text(), parse_constant=lambda _: None)
    return {"generated_from": str(path), **_jsonable(results)}


@app.get("/", include_in_schema=False)
def ui() -> FileResponse:
    """The demo UI: pick a strain, run the funnel, see the ranking and the benchmark.

    A single self-contained file served same-origin, so it talks to the endpoints
    above with no CORS and no build step. It is a thin view over the API — every
    number it shows comes from a documented endpoint, nothing is computed here.
    """
    return FileResponse(STATIC_DIR / "index.html")
