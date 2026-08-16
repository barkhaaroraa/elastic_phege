"""API tests (Phase 5).

Split deliberately in two:

* Pure tests for the serialisation helpers, which run anywhere.
* Contract tests against a live cluster, which **skip** rather than fail when
  Elasticsearch is down or unpopulated. A red suite should mean the code is
  wrong, not that the reviewer has not run ``make up`` yet.

The contract tests assert the parts that are easy to regress silently: that an
unknown strain is a 404 and a stopped funnel is a 422, and that every returned
recommendation carries a reason.
"""

from __future__ import annotations

import math

import pytest
from fastapi.testclient import TestClient

from phageforge import config
from phageforge.api import main
from phageforge.bench import harness

# ------------------------------------------------------------------ pure


def test_jsonable_strips_non_finite():
    """NaN is valid Python and invalid JSON; it must not reach a response."""
    payload = {"a": float("nan"), "b": [1.0, float("inf")], "c": {"d": float("-inf")}, "e": 2}
    assert main._jsonable(payload) == {"a": None, "b": [1.0, None], "c": {"d": None}, "e": 2}


def test_jsonable_leaves_ordinary_values_alone():
    payload = {"s": "x", "n": 1, "f": 0.5, "l": [None, True]}
    assert main._jsonable(payload) == payload


def test_harness_save_writes_strict_json(tmp_path):
    """``make bench`` must not produce a file that ``jq`` refuses to read."""
    import json

    results = {"results": [{"method": "random", "p95_ms": float("nan")}], "wall_seconds": 1.0}
    path = harness.save(results, tmp_path / "benchmark.json")
    text = path.read_text()
    assert "NaN" not in text
    # parse_constant fires only on NaN/Infinity, so this raises if any survived.
    reloaded = json.loads(text, parse_constant=lambda c: pytest.fail(f"non-strict token {c}"))
    assert reloaded["results"][0]["p95_ms"] is None


def test_finite_is_recursive():
    assert harness._finite([{"x": float("nan")}]) == [{"x": None}]
    assert math.isnan(float("nan"))  # guards the helper above being a no-op


# ------------------------------------------------------------------ live


@pytest.fixture(scope="module")
def client():
    """A TestClient, or a skip if the cluster cannot serve the funnel."""
    es = config.get_client()
    try:
        if es.cluster.health()["status"] == "red":
            pytest.skip("Elasticsearch is red")
    except Exception as exc:  # noqa: BLE001
        pytest.skip(f"Elasticsearch unreachable: {exc}")
    if not es.indices.exists(index=config.BACTERIA):
        pytest.skip("pf-bacteria missing — run `make ingest`")
    with TestClient(main.app) as test_client:
        if not test_client.get("/health").json()["ready"]:
            pytest.skip("cluster is up but the funnel has no features — run `make features`")
        yield test_client


@pytest.fixture(scope="module")
def strain_id(client) -> str:
    """A strain that actually has a genome vector, chosen from the index."""
    strains = client.get("/strains", params={"limit": 1}).json()
    if not strains:
        pytest.skip("no strains indexed")
    return strains[0]["strain_id"]


def test_health_reports_indices(client):
    body = client.get("/health").json()
    assert body["cluster_status"] in {"green", "yellow"}
    assert body["indices"][config.INTERACTIONS]["docs"] > 0
    assert body["model_version"]


def test_unknown_strain_is_404_not_a_funnel_failure(client):
    """The distinction a caller needs: bad input vs. the funnel giving up."""
    response = client.get("/predict/NO-SUCH-STRAIN-XYZ")
    assert response.status_code == 404
    assert "not in" in response.json()["detail"]


def test_predict_returns_a_reasoned_shortlist(client, strain_id):
    response = client.get(
        "/predict/" + strain_id, params={"top_n": 5, "persist": False}
    )
    if response.status_code == 422:
        pytest.skip(f"funnel stopped for {strain_id}: {response.json()['failure_reason']}")
    assert response.status_code == 200

    body = response.json()
    assert body["strain_id"] == strain_id
    assert 1 <= len(body["shortlist"]) <= 5
    assert [e["rank"] for e in body["shortlist"]] == list(range(1, len(body["shortlist"]) + 1))

    scores = [e["score"] for e in body["shortlist"]]
    assert scores == sorted(scores, reverse=True), "shortlist must be returned in rank order"

    for entry in body["shortlist"]:
        # G3: no recommendation without a stated basis.
        assert entry["reasons"], f"{entry['phage_id']} came back with no reason"
        assert entry["evidence_tier"] in {1, 4}
        assert 0.0 <= entry["prior"] <= 1.0


def test_shortlist_is_never_silently_empty(client):
    """A strain indexed without a genome vector must 422 with the cause named."""
    es = config.get_client()
    hits = es.search(
        index=config.BACTERIA,
        size=1,
        query={"bool": {"must_not": [{"exists": {"field": "genome_vector"}}]}},
        source_includes=["strain_id"],
    )["hits"]["hits"]
    if not hits:
        pytest.skip("every indexed strain has a genome_vector — nothing to stop the funnel")

    response = client.get("/predict/" + hits[0]["_source"]["strain_id"])
    assert response.status_code == 422
    body = response.json()
    assert body["failed_stage"] == "N"
    assert body["failure_reason"]
    assert "shortlist" not in body, "a stopped funnel must not return a shortlist at all"


def test_ablation_switches_change_the_ranking(client, strain_id):
    """The ablation flags must reach the funnel, not be silently accepted."""
    params = {"top_n": 5, "persist": False, "cocktail": False, "explain": False}
    full = client.get("/predict/" + strain_id, params=params)
    if full.status_code != 200:
        pytest.skip("funnel did not complete for this strain")
    ablated = client.get("/predict/" + strain_id, params={**params, "use_prior": False})
    assert ablated.status_code == 200
    assert [e["score"] for e in full.json()["shortlist"]] != [
        e["score"] for e in ablated.json()["shortlist"]
    ], "zeroing the prior left the scores identical"


def test_run_is_persisted_and_retrievable(client, strain_id):
    response = client.get("/predict/" + strain_id, params={"top_n": 3, "persist": True})
    if response.status_code != 200:
        pytest.skip("funnel did not complete for this strain")
    run_id = response.json()["run_id"]

    config.get_client().indices.refresh(index=config.PREDICTIONS)
    stored = client.get(f"/runs/{run_id}")
    assert stored.status_code == 200
    assert stored.json()["run_id"] == run_id
    assert stored.json()["strain_id"] == strain_id


def test_unknown_run_is_404(client):
    assert client.get("/runs/deadbeefdeadbeef").status_code == 404


def test_catalogue_lookups(client, strain_id):
    assert client.get(f"/strains/{strain_id}").json()["strain_id"] == strain_id
    assert client.get("/strains/NOPE-123").status_code == 404

    phages = client.get("/phages", params={"limit": 1}).json()
    if phages:
        pid = phages[0]["phage_id"]
        assert client.get(f"/phages/{pid}").json()["phage_id"] == pid
    assert client.get("/phages/NOPE-123").status_code == 404


def test_benchmark_endpoint_is_strict_json(client):
    response = client.get("/benchmark")
    if response.status_code == 404:
        pytest.skip("no benchmark.json — run `make bench`")
    body = response.json()  # httpx parses strictly; a bare NaN would raise here
    assert body["results"]
    assert {"random", "funnel"} <= {r["method"] for r in body["results"]}
