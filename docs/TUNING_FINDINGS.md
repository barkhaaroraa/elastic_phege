# PhageForge — Ranking Accuracy & Latency Findings

**Created:** 2026-08-16 · **Scope:** measured answers to "make it more accurate and
cheaper", on the current index

Everything here is measured on the same 390 strains, the same 5 strain-grouped
folds and the same seed (20260815) the published benchmark uses. Nothing is
argued from first principles.

> **Status:** findings only. No production file was edited — `funnel/stages.py`
> and `bench/harness.py` were being modified by a parallel session while this ran
> (§6). The changes below are proposals with the numbers already attached.

---

## 1. Headline

| | P@10 | R@10 | AUPRC | median ms | p95 ms |
|---|--:|--:|--:|--:|--:|
| funnel as shipped | 0.6333 | 0.4130 | 0.6548 | 58.9 | 65.2 |
| **with the four changes in §3** | **0.6749** | **0.4411** | **0.6951** | **25.8** | **30.2** |
| delta | **+0.0415** | +0.0281 | **+0.0403** | **−56%** | **−54%** |

Paired bootstrap over the 390 strains, 10,000 resamples:

- ΔP@10 = **+0.0415**, 95% CI **[+0.0300, +0.0544]** — 136 strains better, 56 worse
- ΔAUPRC = **+0.0403**, 95% CI **[+0.0279, +0.0526]** — 276 better, 110 worse

Hyperparameters were chosen by **nested selection**: for each fold, the
configuration was picked on the other four folds and applied to the held-out one.
All five folds independently selected the same point, so the nested and
full-data numbers coincide.

**Why this matters beyond the metric.** `REVIEW_BACKLOG.md` §2.1 flags that the
funnel beats `phylo_nn` (0.585) by only 4.8 points, and calls that the number a
reviewer will press on. This takes the margin to **+9.0 points** — it roughly
doubles the only result that distinguishes the funnel from a twenty-line
baseline, and it does so while getting faster.

### 1.1 How much room is actually left

P@10 is capped per strain at `min(1, positives/10)`; a strain with 5 measured
positives cannot exceed 0.5. Averaged over the 390 evaluable strains that
ceiling is **0.862**.

| | P@10 | % of achievable |
|---|--:|--:|
| shipped | 0.6333 | 73.5% |
| proposed | 0.6749 | 78.3% |

Worth carrying into any write-up: 13.8% of the apparent "error" is a metric
ceiling, not a model failure. The 57 strains scoring P@10 ≤ 0.2 average 5.6
positives against an overall mean of 20.5 — the model is not failing hardest on
hard strains so much as on *thin* ones.

---

## 2. Method

An offline replica of stages N / B / D was built from the indexed data: exact
Euclidean kNN over the MDS coordinates in place of HNSW, the same Gaussian
kernel, the same weighted-average prior, the same per-fold `breadth`. It
reproduces the shipped funnel to four decimal places —

```
shipped funnel (captured from run_funnel)   P@10 0.6331  R@10 0.4129  AUPRC 0.6548
offline replica (exact kNN)                 P@10 0.6333  R@10 0.4130  AUPRC 0.6548
```

— and independently reproduces the published `generalist` baseline at all three
infection cuts (0.581 / 0.413 / 0.274, matching `current_status.md` §3.1
exactly). The 0.0002 gap is HNSW approximation, consistent with the re-ingest
drift already documented there.

That replica is what every sweep below runs on. Fold safety is preserved
throughout: no feature for a strain is ever computed from data in that strain's
own fold.

---

## 3. The four changes

Deltas are each measured in isolation against the shipped configuration.

### 3.1 🔴 Score receptor compatibility — **+1.5 pts P@10**

The largest single win, and it confirms the domain reviewer's §3.1 critique
quantitatively. The parallel session has now landed the `receptor_compat`
feature (smoothed infection-rate lift per host attribute); what was missing was
the weight.

| receptor_compat weight | P@10 | AUPRC |
|---|--:|--:|
| 0 (shipped) | 0.6333 | 0.6548 |
| 3.0 *(currently proposed in `harness.RECEPTOR_WEIGHTS`)* | 0.6462 | 0.6723 |
| **5.0** | **0.6487** | **0.6765** |
| **8.0** | **0.6521** | 0.6769 |
| 12.0 | 0.6500 | 0.6716 |
| 25.0 | 0.6274 | 0.6341 |

The variance-matched 3.0 is a defensible a priori choice and it works — it is
just **under-weighted by roughly half**. Optimum is broad and flat between 5 and
12, so this is not a knife-edge. In the full stack 5.0 wins.

Smoothing also matters and is currently set low:

| `PROFILE_SMOOTHING` | P@10 (at weight 12) |
|---|--:|
| 1 | 0.6441 |
| 5 *(current)* | 0.6500 |
| 10 | 0.6551 |
| 20 | 0.6592 |

Breaking the signal down by attribute (using an independent log-enrichment
formulation, so the conclusion does not depend on one feature definition):
`lps_type` and `capsule_types` carry it, `o_antigen` adds a little, and
**`st` / `phylogroup` / `h_type` actively hurt** — they are proxies for the
phylogeny the prior already encodes, so adding them double-counts it. Keep the
receptor feature to the three receptor fields.

### 3.2 🟠 Shrink the neighbour prior — **+0.5 pts**

The prior is a raw similarity-weighted average. A phage whose evidence rests on
little effective neighbour weight gets the same confidence as one supported by
the whole neighbourhood. Shrinking toward the phage's own fold breadth,
`prior = (Σwy + α·breadth) / (Σw + α)`:

| α | 0 | 0.5 | 1 | 2 | **4** | 8 |
|---|--:|--:|--:|--:|--:|--:|
| P@10 | 0.6333 | 0.6338 | 0.6356 | 0.6364 | **0.6387** | 0.6374 |

This is the fix for the precision-by-confidence curve: top-10 entries at prior
0.2–0.4 are right 40% of the time, at prior >0.8 they are right 86%. The prior
is well calibrated but not sharp, and shrinkage sharpens it.

### 3.3 🟠 Narrow the kernel bandwidth — **+0.3 pts**

`_kernel_weights` sets σ to the median neighbour distance. Multiplying that by
0.75 is better at both metrics:

| σ × | 0.25 | 0.5 | **0.75** | 1.0 *(current)* | 1.5 | 4.0 |
|---|--:|--:|--:|--:|--:|--:|
| P@10 | 0.6213 | 0.6351 | **0.6362** | 0.6333 | 0.6303 | 0.6272 |
| AUPRC | 0.6374 | 0.6569 | **0.6606** | 0.6548 | 0.6448 | 0.6363 |

Neighbour count `k` was swept too (5 → 150). k=15 is marginally best (0.6364)
but the whole range 10–35 is within noise of k=25. **Leave `NEIGHBOUR_K` alone** —
there is no real gain there and changing it invalidates a published number for
nothing.

### 3.4 🔴 Drop Stage A from the score — **+0.4 pts and −33 ms**

Stage A's only route into the ranking is `candidate_rank`. Removing that feature
**improves** accuracy:

| | P@10 | AUPRC |
|---|--:|--:|
| shipped (with `candidate_rank`) | 0.6333 | 0.6548 |
| `candidate_rank` removed | **0.6372** | **0.6613** |

On the improved configuration it is neutral-to-positive (0.6669 → 0.6672). So at
96 phages Stage A is not merely inert as `current_status.md` §5.3 says — it costs
**49% of end-to-end latency** and slightly hurts the ranking. See §4.

This does not argue for deleting Stage A. It argues for **gating it on corpus
size**: skip it when the corpus fits inside `CANDIDATE_POOL`, which is exactly
the condition under which it cannot narrow anything. Stage A's justification is
the 100K–1M tiers (WS-3), and that case is untouched by this.

---

## 4. Efficiency

12 Elasticsearch round trips per prediction, ~90 ms of ES time with `explain`
and `cocktail` on. Measured over 60 strains:

| call | per run | mean ms | ms/prediction |
|---|--:|--:|--:|
| `pf-interactions` agg (Stage B prior, Stage C hosts) | 2.0 | 13.05 | 26.10 |
| `pf-phages` kNN (`_knn_scores`, per-arm) | 2.0 | 7.24 | 14.48 |
| `pf-phages` fetch (Stage A seeds, cocktail vectors) | 2.0 | 6.17 | 12.34 |
| `pf-phages` rrf (Stage A) | 1.0 | 11.54 | 11.54 |
| `msearch` (Stage C significant_terms) | 1.0 | 8.13 | 8.13 |
| `pf-phages` script_score (Stage D) | 1.0 | 6.40 | 6.40 |
| `pf-bacteria` kNN (Stage N) | 1.0 | 5.60 | 5.60 |
| `pf-bacteria` fetch (`get_strain`) | 1.0 | 3.26 | 3.26 |
| `count(pf-phages)` | 1.0 | 1.96 | 1.96 |
| **total** | **12.0** | | **89.8** |

By stage, in benchmark mode (58.9 ms median): N 9.4 ms · B 8.4 ms · **A 32.2 ms
(49%)** · D 8.6 ms · 6.9 ms orchestration.

Measured, cumulatively, by disabling each in turn:

| | median | p95 |
|---|--:|--:|
| as shipped | 58.9 ms | 65.2 ms |
| − the two per-arm kNN searches | 46.1 ms | 55.9 ms |
| − `count(pf-phages)` | 46.2 ms | 53.9 ms |
| − Stage A entirely | **25.8 ms** | **30.2 ms** |

Three things to take from this:

1. **`count(pf-phages)` is a round trip for a display string.** It fills
   `corpus_size` in Stage A's metadata and nothing reads it for scoring. It is a
   constant per index; cache it. 2 ms for free.
2. **Round-trip overhead dominates.** A `count` against a 96-document index costs
   1.96 ms, which is essentially pure transport. Twelve trips is ~25 ms of the
   ~90 ms before any work happens. Stage N + Stage B are strictly sequential
   (B needs N's weights), but `get_strain`, the Stage A seed fetch and the
   cocktail's vector fetch are all independent and could be batched into one
   `msearch`.
3. **The cocktail re-fetches phage documents Stage D already fetched.** Passing
   `rbp_match_vector` through from Stage D's `source_includes` removes a whole
   round trip (~6 ms).

`pf-phages` also carries 89 deleted documents against 96 live ones; a
`_forcemerge` after ingest is worth doing, though at this size it is noise.

---

## 5. Three corrections to claims currently in the docs

### 5.1 🔴 "The RBP arm contributes nothing" — right conclusion, wrong reason

`current_status.md` §3 reading 3 and `REVIEW_BACKLOG.md` §2.2 both attribute the
null RBP ablation to there being nothing to narrow at 96 phages. That is not why.
In the code as benchmarked, `rbp_similarity` and `genome_similarity` were
computed by Stage A, declared in `FEATURE_NAMES`, weighted in `DEFAULT_WEIGHTS`
— and **never referenced by the Painless script**, and never placed in the
`features` map passed to it. The ablation could only ever have measured Stage A's
effect on `candidate_rank`.

The parallel session has already found and fixed this (`_score_script()`
generates the script from `FEATURE_NAMES`, with a test to keep them from
drifting). What the fix does not yet say is what the feature is worth. Wiring it
in and giving it weight makes things **worse**, monotonically:

| `rbp_similarity` weight | 0.0 | 0.5 | 1.0 | 2.0 |
|---|--:|--:|--:|--:|
| P@10 | 0.6749 | 0.6608 | 0.6521 | 0.6290 |

`genome_similarity` is flat (0.6749 → 0.6677 at weight 2.0). So the honest
statement is stronger than the current one and worth stating plainly: the ESM-2
RBP centroid similarity, as computed, is **not** uninformative-because-untested
— it is anti-correlated with infection at this scale. Keeping both weights at
0.0 is correct, and now it is correct for a measured reason.

### 5.2 🟠 `defence_compat` at weight 7.0 is not supported

The variance-matched weights give defence a *larger* coefficient than receptor
(7.0 vs 3.0), because its spread is smaller. But the signal is not there:

| `defence_compat` weight | 0 | 2 | 5 | **7** | 12 | 20 |
|---|--:|--:|--:|--:|--:|--:|
| P@10 | 0.6333 | 0.6344 | 0.6344 | **0.6338** | 0.6287 | 0.6115 |
| AUPRC | 0.6548 | 0.6562 | 0.6543 | **0.6538** | 0.6436 | 0.6213 |

Null at small weight, harmful at large. An independent log-enrichment
formulation reproduces this. Combined with receptor it *reduces* the receptor-only
result (0.6462 → 0.6459 P@10, AUPRC 0.6723 → 0.6673).

This is not a surprise and it should be reported as a **confirmed negative**, not
buried: `REVIEW_BACKLOG.md` §3.2 already states the reason — DefenseFinder gives
*system presence*, and CRISPR/RM presence is not resistance to a *particular*
phage without spacer or motif matching. This measurement is direct evidence for
that caveat, and it strengthens the case for WS-2 (real spacer matching) rather
than weakening it. Recommend shipping `defence_compat` at **0.0** with this table
as the justification.

### 5.3 🟡 `support` is a dead feature

`support = log1p(n_tested)` carries weight 0.30. But the interaction matrix is
dense: across every strain and candidate, `n_tested` only ever takes the values
22–25. Removing the feature entirely is **bit-identical** (0.6333 → 0.6333).

It costs nothing and it is not wrong — but it is one of six declared features and
it can never move a ranking on this dataset. Anything trained on this feature
vector (WS-4) will fit a meaningless coefficient to it. Worth either dropping or
documenting as structurally inert here.

---

## 6. Robustness

### 6.1 Across the infection threshold

Re-derived from the indexed raw 0–4 lysis scores, so no re-ingest was needed and
the vectors were not disturbed. The improvement holds at every cut:

| Cut | base rate | generalist | current scoring | improved (receptor weight tuned per cut) |
|---|--:|--:|--:|--:|
| `>0` | 0.208 | 0.5810 | 0.6372 | **0.6703** |
| `>1` (≥2) | 0.129 | 0.4133 | 0.4845 | **0.5128** |
| `>2` (≥3) | 0.069 | 0.2737 | 0.3448 | **0.3684** |

All P@10. The `generalist` column reproduces the published 0.581 / 0.413 / 0.274
exactly, which is the validation that this offline path is faithful.

**One caveat that must travel with this:** the receptor weight has to come down
as the cut tightens — optimum is ~1.6 at `>0`, ~0.8 at `>1`, ~0.4 at `>2` (in
log-enrichment units). At a fixed weight tuned for `>0`, AUPRC *regresses* at
`>2` (0.4628 → 0.4478). Tuned per cut, both metrics improve at all three. The
mechanism is obvious — the enrichment is estimated from fewer positives — and it
is the right argument for making this weight fitted rather than fixed (WS-4).

### 6.2 What did not work

Recorded so nobody spends the time twice.

- **`defence_compat`** — §5.2.
- **`st` / `phylogroup` / `h_type` as compatibility features** — 0.6323 vs 0.6456
  for the three receptor fields alone. They re-encode the phylogeny.
- **Latent-factor prior (truncated SVD of the training matrix, target factor
  estimated from neighbours)** — works, but only +0.5 pts alone (0.6364 at rank
  32), and it is redundant with the receptor feature: added on top of the full
  stack it *costs* 0.9 pts (0.6669 → 0.6577). Not worth the machinery.
- **Neighbour count `k`** — flat across 10–35. §3.3.
- **`rbp_similarity`** — §5.1.

---

## 7. Recommended sequence

1. **Raise `receptor_compat` to 5.0 and `PROFILE_SMOOTHING` toward 10–20**;
   ship `defence_compat` at 0.0 with §5.2 as the reason. *(Owner: whoever owns
   `stages.py` / `harness.py` right now.)*
2. **Shrink the prior (α = 4) and narrow the kernel (σ × 0.75).** Both are local
   changes to `stage_b_prior` / `_kernel_weights`.
3. **Gate Stage A on `corpus_size <= CANDIDATE_POOL`.** Biggest latency win by
   far, and it improves the ranking. Keep the stage; skip it when it provably
   cannot narrow.
4. **Cache `count(pf-phages)`; pass the cocktail's vectors through from Stage D.**
5. Then WS-4 (LTR) fits these weights properly instead of by sweep — §6.1 shows
   the receptor weight genuinely wants to be data-dependent.

Steps 1–4 are independent of each other and none of them touches the fold
structure, so each can land and be measured on its own.

### Caveats

- Every number here is on **this index, this organism, 96 phages**. The receptor
  finding in particular is the kind that should be re-checked on Klebsiella
  (WS-7) before it is presented as general.
- Latency was measured on the local single-node cluster with a 1 GB heap. The
  *ratios* should hold; the absolute milliseconds will not.
- §3 deltas are measured in isolation and are not additive — the four together
  give +4.2 pts, more than their sum, mostly because dropping `candidate_rank`
  clears room for the receptor feature.
- Reproduction scripts are in the session scratchpad, not in the repo. If these
  findings are acted on, the sweeps belong in `bench/` as a committed experiment
  rather than as a transcript.
