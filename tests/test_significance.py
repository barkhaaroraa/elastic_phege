"""Paired significance tests. Pure — no cluster, no fixtures.

The properties asserted here are the ones that make a reported interval
trustworthy: that pairing is by strain and not by position, that a zero
difference reads as "not detectable" rather than "significant", and that a real
difference is actually detected at the sample size the benchmark runs at.
"""

from __future__ import annotations

import math

import pytest

from phageforge.bench import significance


def _scores(values: list[float], prefix: str = "s") -> dict[str, float]:
    return {f"{prefix}{i}": v for i, v in enumerate(values)}


# ------------------------------------------------------------------ align


def test_align_pairs_on_strain_id_not_position():
    """The whole point of recording strain ids.

    ``b`` is missing s1, so a positional zip would pair s2's score against s3's
    and report a difference that never happened.
    """
    a = {"s1": 1.0, "s2": 2.0, "s3": 3.0}
    b = {"s2": 2.0, "s3": 3.0}
    ids, x, y = significance.align(a, b)
    assert ids == ["s2", "s3"]
    assert list(x) == [2.0, 3.0]
    assert list(y) == [2.0, 3.0]


def test_align_drops_nan_pairs():
    """Recall is NaN on a strain with no positives; it cannot enter a difference."""
    a = {"s1": 0.5, "s2": float("nan")}
    b = {"s1": 0.4, "s2": 0.9}
    ids, x, y = significance.align(a, b)
    assert ids == ["s1"]
    assert len(x) == len(y) == 1


def test_align_of_disjoint_methods_is_empty():
    ids, x, y = significance.align({"a": 1.0}, {"b": 1.0})
    assert ids == []
    assert len(x) == len(y) == 0


# ---------------------------------------------------------------- compare


def test_identical_methods_are_not_significant():
    """Two configurations differing only in a zero-weighted feature.

    This is the exact case the wiring fix had to land safely: every strain tied,
    so the honest verdict is "no detectable difference", not p < 0.05.
    """
    scores = _scores([0.1, 0.4, 0.7, 0.2, 0.9, 0.3])
    result = significance.compare(scores, dict(scores), method="a", reference="b")
    assert result is not None
    assert result.difference == 0.0
    assert result.p_wilcoxon == 1.0
    assert not result.significant
    assert result.n_tied == 6
    assert "no detectable difference" in result.verdict()


def test_a_tiny_gap_at_realistic_n_is_not_significant():
    """The 0.0002 that retired the RBP arm, reproduced.

    390 strains, one of which differs by a single P@10 step. A mean would report
    a non-zero gap; the interval must straddle zero.
    """
    base = [0.6] * 390
    nudged = list(base)
    nudged[0] = 0.7
    result = significance.compare(
        _scores(nudged), _scores(base), method="funnel", reference="funnel_minus_rbp"
    )
    assert result is not None
    assert result.difference == pytest.approx(0.1 / 390)
    assert not result.significant, "a one-strain difference must not read as real"


def test_a_real_gap_at_realistic_n_is_significant():
    """+0.05 on most strains at n=390 must be detected, or the test is useless."""
    reference = [0.55] * 390
    method = [0.60] * 390
    result = significance.compare(
        _scores(method), _scores(reference), method="funnel", reference="phylo_nn"
    )
    assert result is not None
    assert result.significant
    assert result.ci_low > 0
    assert result.difference == pytest.approx(0.05)


def test_confidence_interval_brackets_the_observed_difference():
    reference = [0.3, 0.5, 0.4, 0.6, 0.2, 0.7, 0.5, 0.4]
    method = [0.5, 0.6, 0.5, 0.8, 0.3, 0.7, 0.6, 0.5]
    result = significance.compare(_scores(method), _scores(reference), method="a", reference="b")
    assert result is not None
    assert result.ci_low <= result.difference <= result.ci_high


def test_difference_sign_follows_the_method_not_the_reference():
    """A method worse than its reference must report a negative difference."""
    result = significance.compare(
        _scores([0.2] * 50), _scores([0.6] * 50), method="worse", reference="better"
    )
    assert result is not None
    assert result.difference < 0
    assert result.ci_high < 0
    assert "worse" in result.verdict()


def test_win_loss_tie_counts_sum_to_the_pair_count():
    method = _scores([0.5, 0.6, 0.4, 0.5])
    reference = _scores([0.4, 0.6, 0.6, 0.5])
    result = significance.compare(method, reference, method="a", reference="b")
    assert result is not None
    assert result.n_better + result.n_worse + result.n_tied == result.n_pairs == 4
    assert (result.n_better, result.n_worse, result.n_tied) == (1, 1, 2)


def test_comparison_is_reproducible_across_calls():
    """A seeded bootstrap. An unseeded interval would drift between reruns."""
    method, reference = _scores([0.3, 0.5, 0.7, 0.2, 0.9]), _scores([0.2, 0.5, 0.6, 0.3, 0.7])
    first = significance.compare(method, reference, method="a", reference="b")
    second = significance.compare(method, reference, method="a", reference="b")
    assert first is not None and second is not None
    assert (first.ci_low, first.ci_high) == (second.ci_low, second.ci_high)


def test_too_few_pairs_returns_none_rather_than_a_fake_interval():
    assert significance.compare({"s1": 0.5}, {"s1": 0.4}, method="a", reference="b") is None


# ----------------------------------------------------------- compare_all


def test_compare_all_excludes_the_reference_and_sorts_by_difference():
    per_strain = {
        "phylo_nn": _scores([0.5] * 40),
        "funnel": _scores([0.6] * 40),
        "random": _scores([0.2] * 40),
    }
    report = significance.compare_all(per_strain, reference="phylo_nn")
    assert [c.method for c in report.comparisons] == ["funnel", "random"]
    assert report.comparisons[0].difference > report.comparisons[1].difference


def test_compare_all_rejects_an_unknown_reference():
    with pytest.raises(KeyError, match="nonexistent"):
        significance.compare_all({"a": _scores([0.1] * 5)}, reference="nonexistent")


def test_report_round_trips_through_its_dict_form():
    """``format_table`` rebuilds Comparisons from JSON; the keys must line up."""
    per_strain = {"funnel": _scores([0.6] * 30), "phylo_nn": _scores([0.5] * 30)}
    report = significance.compare_all(per_strain, reference="phylo_nn")
    payload = report.to_dict()
    rebuilt = [
        significance.Comparison(
            **{k: v for k, v in c.items() if k not in ("significant", "verdict")}
        )
        for c in payload["comparisons"]
    ]
    assert [c.difference for c in rebuilt] == [c.difference for c in report.comparisons]


def test_formatted_report_states_what_a_spanning_interval_means():
    per_strain = {"funnel": _scores([0.6] * 30), "phylo_nn": _scores([0.5] * 30)}
    text = significance.format_report(significance.compare_all(per_strain, reference="phylo_nn"))
    assert "not distinguishable from none" in text
    assert "funnel" in text


def test_all_values_are_json_safe():
    """``harness.save`` writes with ``allow_nan=False``; a NaN here would break it."""
    per_strain = {"a": _scores([0.5] * 20), "b": _scores([0.5] * 20)}
    payload = significance.compare_all(per_strain, reference="b").to_dict()
    for comparison in payload["comparisons"]:
        for key, value in comparison.items():
            if isinstance(value, float):
                assert math.isfinite(value), f"{key} is not finite"
