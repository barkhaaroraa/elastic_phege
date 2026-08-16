"""Stage D scoring contract, and the pure helpers the stages are built on.

Split like ``test_api.py``: pure tests that run anywhere, and contract tests that
**skip** when the cluster is down or unpopulated.

The centre of gravity here is the agreement between three things that must all
name the same features:

* ``FEATURE_NAMES``      -- the declared feature vector
* ``DEFAULT_WEIGHTS``    -- a weight per declared feature
* the ``features`` dict  -- what ``stage_d_rank`` actually builds per candidate

They drifted once, silently, and it cost a wrong conclusion: ``rbp_similarity``
and ``genome_similarity`` were declared and weighted but never referenced by the
Painless source, so the "- RBP arm" ablation could only act through
``candidate_rank`` and reported that the RBP arm changed nothing. That was a fact
about the wiring being read as a fact about the biology. These tests exist so the
same drift fails loudly instead.
"""

from __future__ import annotations

import math

import pytest

from phageforge import config
from phageforge.funnel import stages

# ------------------------------------------------------------------ pure


def test_every_declared_feature_has_a_weight():
    missing = set(stages.FEATURE_NAMES) - set(stages.DEFAULT_WEIGHTS)
    assert not missing, f"declared in FEATURE_NAMES but unweighted: {sorted(missing)}"


def test_no_weight_without_a_declared_feature():
    """``bias`` is the one legitimate weight with no matching feature."""
    extra = set(stages.DEFAULT_WEIGHTS) - set(stages.FEATURE_NAMES) - {"bias"}
    assert not extra, f"weighted but never scored: {sorted(extra)}"


def test_score_script_references_every_declared_feature():
    """The generated Painless source must consume the whole feature vector.

    This is the assertion that would have caught the original bug on the day it
    was introduced.
    """
    script = stages._score_script()
    for name in stages.FEATURE_NAMES:
        assert f"params.w.{name}" in script, f"{name} is declared but never scored"
        assert f"f.{name}" in script, f"{name}'s value is never read"


def test_score_script_is_generated_not_hardcoded():
    """Adding a feature must change the script without anyone editing Painless."""
    original = stages.FEATURE_NAMES
    try:
        stages.FEATURE_NAMES = (*original, "a_new_feature")
        assert "params.w.a_new_feature" in stages._score_script()
    finally:
        stages.FEATURE_NAMES = original


def test_cosine_features_are_declared_features():
    unknown = stages._COSINE_FEATURES - set(stages.FEATURE_NAMES)
    assert not unknown, f"centred but not scored: {sorted(unknown)}"


def test_distance_from_score_inverts_the_l2_norm_formula():
    """ES scores ``l2_norm`` as ``1 / (1 + d^2)``; the inverse must round-trip."""
    for distance in (0.0, 1e-4, 0.05, 0.117, 2.0):
        score = 1.0 / (1.0 + distance**2)
        assert stages._distance_from_score(score) == pytest.approx(distance, abs=1e-9)


def test_distance_from_score_handles_a_zero_score():
    assert stages._distance_from_score(0.0) == float("inf")


def test_kernel_weights_are_monotonic_in_distance():
    """A nearer strain must never be weighted below a further one."""
    weights = stages._kernel_weights([0.0, 0.01, 0.05, 0.1, 0.2])
    assert weights == sorted(weights, reverse=True)
    assert weights[0] == pytest.approx(1.0)


def test_kernel_bandwidth_adapts_to_the_neighbourhood():
    """Scaling every distance must leave the weights unchanged.

    ``sigma`` is the median neighbour distance, so the kernel is scale-free: a
    tight ST cluster and a sparse corner of the phylogeny get the same shape.
    A fixed bandwidth would collapse one of the two.
    """
    tight = stages._kernel_weights([0.001, 0.002, 0.004])
    loose = stages._kernel_weights([0.1, 0.2, 0.4])
    assert tight == pytest.approx(loose)


def test_kernel_weights_survive_an_all_zero_neighbourhood():
    """Duplicate genomes give distance 0 throughout; that must not divide by zero."""
    assert stages._kernel_weights([0.0, 0.0]) == [1.0, 1.0]


def test_weighted_centroid_is_unit_length():
    centroid = stages._weighted_centroid([[3.0, 0.0], [0.0, 4.0]], [1.0, 1.0])
    assert math.isclose(math.sqrt(sum(v * v for v in centroid)), 1.0)


def test_weighted_centroid_of_nothing_is_none():
    assert stages._weighted_centroid([], []) is None


def test_weighted_centroid_of_opposed_vectors_is_none():
    """Cancelling to the origin has no direction, so there is no query vector."""
    assert stages._weighted_centroid([[1.0, 0.0], [-1.0, 0.0]], [1.0, 1.0]) is None


def test_empty_stage_names_the_stage_and_the_cause():
    exc = stages.EmptyStage("B", "no neighbour has a strain-resolution row")
    assert exc.stage == "B"
    assert "strain-resolution" in exc.reason
    assert "stage B" in str(exc)


# -------------------------------------------------------------- contract


@pytest.fixture(scope="module")
def client():
    """A live client, or a skip. A red suite must mean the code is wrong."""
    try:
        es = config.get_client()
        if es.cluster.health()["status"] == "red":
            pytest.skip("Elasticsearch is red")
    except Exception as exc:  # noqa: BLE001 - unreachable is a skip, not a failure
        pytest.skip(f"Elasticsearch unreachable: {exc}")
    if not es.indices.exists(index=config.BACTERIA):
        pytest.skip("pf-bacteria missing — run `make ingest`")
    if not es.count(index=config.BACTERIA, query={"exists": {"field": "genome_vector"}})["count"]:
        pytest.skip("no genome vectors — run `make features`")
    return es


@pytest.fixture(scope="module")
def scored(client):
    """Run N -> B -> A -> D once and hand back the candidates and the ranking."""
    strain_id = client.search(
        index=config.BACTERIA,
        size=1,
        query={"exists": {"field": "genome_vector"}},
        source_includes=["strain_id"],
    )["hits"]["hits"][0]["_source"]["strain_id"]

    strain = stages.get_strain(client, strain_id)
    neighbours = stages.find_neighbours(client, strain)
    priors = stages.stage_b_prior(client, neighbours.results)
    try:
        candidates = stages.stage_a_candidates(client, priors.results)
    except stages.EmptyStage as exc:
        pytest.skip(f"stage A produced nothing for {strain_id}: {exc.reason}")
    ranked = stages.stage_d_rank(
        client, candidates.results, priors.results, top_n=96, strain=strain
    )
    return candidates, ranked


def test_stage_d_builds_exactly_the_declared_features(scored):
    """The drift guard, against a real run rather than the source.

    A feature can be declared, weighted, and emitted into the Painless source and
    *still* never reach the score if ``stage_d_rank`` does not build it -- the
    script's ``containsKey`` guard would quietly substitute 0. This is the test
    that catches that.
    """
    _, ranked = scored
    for entry in ranked.results:
        built = set(entry["features"])
        assert built == set(stages.FEATURE_NAMES), (
            f"{entry['phage_id']}: built {sorted(built)}, "
            f"declared {sorted(stages.FEATURE_NAMES)}"
        )


def test_stage_a_similarities_are_never_below_orthogonal(scored):
    """ES reports cosine as ``(1 + cos) / 2``; an unretrieved candidate is 0.5.

    A default of 0.0 would mean *anti-similar* and rank an unretrieved phage below
    every retrieved one -- invisible at 96 phages where the pool exceeds the
    corpus, a silent misranking at the tiers Stage A exists for.
    """
    candidates, _ = scored
    for candidate in candidates.results:
        for arm in ("rbp_similarity", "genome_similarity"):
            assert 0.0 <= candidate[arm] <= 1.0
            if candidate[arm] == 0.0:
                pytest.fail(f"{candidate['phage_id']} scored {arm} at the anti-similar floor")


def test_cosine_features_reach_stage_d_centred(scored):
    """Centring on 0.5 makes a zero weight genuinely neutral."""
    candidates, ranked = scored
    raw = {c["phage_id"]: c for c in candidates.results}
    for entry in ranked.results:
        for name in stages._COSINE_FEATURES:
            assert entry["features"][name] == pytest.approx(raw[entry["phage_id"]][name] - 0.5)


def test_zero_weighted_features_do_not_move_the_score(scored):
    """The invariant that let the wiring fix land on its own.

    Every feature added with weight 0.0 must leave the ranking bit-identical, so a
    plumbing change and a scoring change are never entangled in one commit.
    """
    _, ranked = scored
    zeroed = [n for n in stages.FEATURE_NAMES if stages.DEFAULT_WEIGHTS.get(n) == 0.0]
    if not zeroed:
        pytest.skip("no zero-weighted features to check")
    for entry in ranked.results:
        expected = stages.DEFAULT_WEIGHTS["bias"] + sum(
            stages.DEFAULT_WEIGHTS[name] * entry["features"][name]
            for name in stages.FEATURE_NAMES
            if name not in zeroed
        )
        assert entry["score"] == pytest.approx(max(expected, 0.0), abs=1e-4)
