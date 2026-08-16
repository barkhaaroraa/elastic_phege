"""Paired significance testing for the benchmark.

``harness.py`` reports a mean per method. A mean on its own cannot answer the only
question that matters when two methods are close: **is the gap real?**

The benchmark's own history is the argument for this module. The "- RBP arm"
ablation was reported as `0.6333 vs 0.6331` and read as "the RBP arm contributes
exactly nothing" -- a four-decimal comparison used to retire a whole phase of
work. Re-running the same configuration after an unrelated re-ingest moved the
same pair to `0.6341 vs 0.6331`. Neither gap was ever shown to be distinguishable
from zero, and neither was shown to be distinguishable from the other.

Every method is scored on the *same* held-out strains, so the comparisons here are
**paired**: the unit of analysis is one strain's difference between two methods,
not two independent means. Pairing removes between-strain variance, which is by
far the largest source of spread -- a strain with 40 susceptible phages scores
better under every method than a strain with 2.

Two tests, deliberately:

``wilcoxon``
    Signed-rank on the paired differences. Distribution-free, which matters
    because per-strain P@10 on 96 phages is discrete (multiples of 0.1) and
    nowhere near normal.
``bootstrap``
    Percentile bootstrap over strains, giving a **confidence interval on the
    difference**. This is the one to quote: "+0.048 P@10 (95% CI +0.021 to
    +0.076)" says everything "+4.8 points" leaves out.

A p-value alone would invite the opposite error -- declaring a 0.001 difference
"significant" at n=390 when it is far too small to matter. Report the interval and
let the reader see both the sign and the size.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

#: Bootstrap resamples. 10,000 puts the Monte-Carlo error on a 95% bound at well
#: under a thousandth of a P@10 point, which is finer than any difference the
#: benchmark can meaningfully report.
N_RESAMPLES = 10_000

#: Fixed so a reported interval is reproducible. The benchmark already learned
#: that unseeded variation gets mistaken for signal.
SEED = 20260815


@dataclass
class Comparison:
    """One method measured against a reference, paired by strain."""

    method: str
    reference: str
    metric: str
    n_pairs: int
    mean_method: float
    mean_reference: float
    difference: float
    ci_low: float
    ci_high: float
    p_wilcoxon: float
    n_better: int
    n_worse: int
    n_tied: int

    @property
    def significant(self) -> bool:
        """A CI excluding zero. Equivalent to the interval not straddling no-effect."""
        return (self.ci_low > 0.0) or (self.ci_high < 0.0)

    def verdict(self) -> str:
        """One line a reader can act on, sign and size included."""
        if not self.significant:
            return (
                f"no detectable difference ({self.difference:+.4f}, "
                f"95% CI {self.ci_low:+.4f} to {self.ci_high:+.4f})"
            )
        direction = "better" if self.difference > 0 else "worse"
        return (
            f"{abs(self.difference):.4f} {direction} "
            f"(95% CI {self.ci_low:+.4f} to {self.ci_high:+.4f}, p={self.p_wilcoxon:.2g})"
        )

    def to_dict(self) -> dict:
        return {
            "method": self.method,
            "reference": self.reference,
            "metric": self.metric,
            "n_pairs": self.n_pairs,
            "mean_method": self.mean_method,
            "mean_reference": self.mean_reference,
            "difference": self.difference,
            "ci_low": self.ci_low,
            "ci_high": self.ci_high,
            "p_wilcoxon": self.p_wilcoxon,
            "significant": self.significant,
            "n_better": self.n_better,
            "n_worse": self.n_worse,
            "n_tied": self.n_tied,
            "verdict": self.verdict(),
        }


def align(
    a: dict[str, float], b: dict[str, float]
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Pair two ``{strain_id: score}`` maps on their shared strains.

    Methods do not always cover the same strains -- a funnel configuration can
    stop on a strain that a baseline ranks fine, and the harness records that as a
    failure and moves on. Comparing the surviving means would then compare two
    different strain populations. Intersecting first is the only honest option,
    and the count is reported so a large drop is visible.

    Strains where either score is NaN (recall on a strain with no positives) are
    dropped for the same reason.
    """
    shared = sorted(set(a) & set(b))
    pairs = [(s, a[s], b[s]) for s in shared]
    usable = [(s, x, y) for s, x, y in pairs if not (math.isnan(x) or math.isnan(y))]
    ids = [s for s, _, _ in usable]
    return ids, np.array([x for _, x, _ in usable]), np.array([y for _, _, y in usable])


def compare(
    method_scores: dict[str, float],
    reference_scores: dict[str, float],
    *,
    method: str,
    reference: str,
    metric: str = "p_at_10",
    n_resamples: int = N_RESAMPLES,
    seed: int = SEED,
) -> Comparison | None:
    """Paired bootstrap CI and Wilcoxon p-value for ``method`` minus ``reference``.

    Returns ``None`` when fewer than two strains are shared, which is not an error
    -- it is the honest answer when there is nothing to compare.
    """
    ids, x, y = align(method_scores, reference_scores)
    if len(ids) < 2:
        return None

    differences = x - y
    rng = np.random.default_rng(seed)
    # Resample *strains*, not differences-and-scores independently: the strain is
    # the sampling unit, so a resample must take a strain's paired difference
    # whole. Drawing indices does exactly that.
    picks = rng.integers(0, len(differences), size=(n_resamples, len(differences)))
    means = differences[picks].mean(axis=1)
    ci_low, ci_high = (float(v) for v in np.percentile(means, [2.5, 97.5]))

    p_value = _wilcoxon(differences)

    return Comparison(
        method=method,
        reference=reference,
        metric=metric,
        n_pairs=len(ids),
        mean_method=float(x.mean()),
        mean_reference=float(y.mean()),
        difference=float(differences.mean()),
        ci_low=ci_low,
        ci_high=ci_high,
        p_wilcoxon=p_value,
        n_better=int((differences > 0).sum()),
        n_worse=int((differences < 0).sum()),
        n_tied=int((differences == 0).sum()),
    )


def _wilcoxon(differences: np.ndarray) -> float:
    """Two-sided signed-rank p-value, or 1.0 when every pair is tied.

    ``scipy.stats.wilcoxon`` raises rather than returning 1.0 when all differences
    are zero. All-tied is a perfectly ordinary outcome here -- two funnel
    configurations that differ only in a zero-weighted feature produce identical
    rankings on every strain -- and it means "no evidence of a difference", which
    is exactly p = 1.
    """
    if not differences.size or not np.any(differences):
        return 1.0
    from scipy.stats import wilcoxon

    return float(wilcoxon(differences, zero_method="wilcox", alternative="two-sided").pvalue)


@dataclass
class Report:
    """Every method compared against one reference, for one metric."""

    reference: str
    metric: str
    comparisons: list[Comparison] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "reference": self.reference,
            "metric": self.metric,
            "comparisons": [c.to_dict() for c in self.comparisons],
        }


def compare_all(
    per_strain: dict[str, dict[str, float]],
    *,
    reference: str,
    metric: str = "p_at_10",
    n_resamples: int = N_RESAMPLES,
    seed: int = SEED,
) -> Report:
    """Compare every method in ``per_strain`` against ``reference``.

    ``per_strain`` is ``{method: {strain_id: score}}`` as written by the harness.
    """
    if reference not in per_strain:
        raise KeyError(f"reference method {reference!r} not in {sorted(per_strain)}")
    report = Report(reference=reference, metric=metric)
    for name in sorted(per_strain):
        if name == reference:
            continue
        comparison = compare(
            per_strain[name],
            per_strain[reference],
            method=name,
            reference=reference,
            metric=metric,
            n_resamples=n_resamples,
            seed=seed,
        )
        if comparison is not None:
            report.comparisons.append(comparison)
    report.comparisons.sort(key=lambda c: -c.difference)
    return report


def format_report(report: Report) -> str:
    """Plain-text table, matching ``harness.format_table``'s house style."""
    lines = [
        f"Paired comparison vs {report.reference} — {report.metric}",
        f"{'Method':30s} {'diff':>8s} {'95% CI':>19s} {'p':>9s} {'W/L/T':>12s}",
        "-" * 82,
    ]
    for c in report.comparisons:
        flag = " *" if c.significant else "  "
        lines.append(
            f"{c.method:30s} {c.difference:+8.4f} "
            f"{c.ci_low:+8.4f}..{c.ci_high:+8.4f} {c.p_wilcoxon:9.2g} "
            f"{c.n_better:>4d}/{c.n_worse:<4d}/{c.n_tied:<4d}{flag}"
        )
    lines.append("")
    lines.append("* 95% CI excludes zero. Paired by strain; W/L/T counts strains.")
    lines.append(
        "A CI spanning zero means the difference is not distinguishable from none — "
        "not that the methods are equal."
    )
    return "\n".join(lines)
