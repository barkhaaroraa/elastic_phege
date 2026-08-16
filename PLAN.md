# Bare metal: close the two functional gaps, then lock the repo down

## Context

PhageForge runs end to end today — ingest → features → proteins → funnel → API → benchmark
— and `data/derived/benchmark.json` carries real strain-grouped numbers (funnel P@10 0.633
vs generalist 0.581, phylo-NN 0.585, random 0.204; p95 74 ms; cocktail coverage 92%).
`current_status.md` is stale and says Phases 4–6 are not started; they are.

Two things are genuinely unfinished, and both are load-bearing for the claims the project
makes about itself:

1. **The RBP ablation measures nothing.** Stage A computes `rbp_similarity` and
   `genome_similarity` per candidate (`stages.py:502`), Stage D drops them on the floor —
   `DEFAULT_WEIGHTS` sets both to `0.0` (`stages.py:694`) and the Painless script
   (`stages.py:772`) never references them. So `funnel_minus_rbp` can only differ through
   Stage A's RRF rank ordering, and the measured gap is 0.0002 *in the ablation's favour*.
   The docstring at `stages.py:479` shows this was noticed once and only half-fixed.
   TECHNICAL_DESIGN §7.6 sells the ablations as "which part carries the signal"; right now
   one of them cannot answer that question.
2. **The hybrid-search half of Phase 6 does not exist.** `current_status.md` §7 lists
   "recall@10 vs exact kNN, p95 latency, RRF-vs-pure-vector" as part of the deliverable, and
   §4.1 makes a testable prediction — that int8/BBQ quantization, which destroyed the
   phylogeny at 404×256, is the right call at scale. Nothing measures either.

The friend's domain review (SOLUTION_OVERVIEW §12) lands after this. Per your call, we are
**not** pre-building config seams for it — current defaults stay hardcoded.

## Ownership boundary

Peer session `0a` is doing a tests-and-polish pass and is actively editing `api/main.py`,
`api/models.py`, `bench/harness.py`, `features/proteins.py`, `funnel/pipeline.py`,
`ingest/gaborieau.py`, `tests/test_api.py`, `pyproject.toml` and `Makefile`.

**This plan touches none of those except two surgical additions.** Files claimed here:

| File | Status |
|---|---|
| `src/phageforge/funnel/stages.py` | exclusive — peer has not touched it |
| `src/phageforge/bench/hybrid.py` | new |
| `src/phageforge/bench/weights.py` | new |
| `tests/test_stages.py`, `tests/test_hybrid.py` | new |
| `README.md` | new |
| `current_status.md` | exclusive — peer is not in docs |
| `src/phageforge/cli.py` | **shared** — add two subparsers only |
| `Makefile` | **shared** — add two targets only |

Do the `cli.py` and `Makefile` edits last, re-reading each file immediately before editing.

---

## Step 0 — git baseline (do this first)

There is no repository at all. Three sessions are editing the same tree with no way to
diff, revert, or review — and the friend's feedback is about to arrive on top.

```bash
git init && git add -A && git commit -m "Baseline: working funnel, API and benchmark"
```

`.gitignore` already excludes `.venv/`, `data/`, caches. Add `.pytest_cache/` and
`.ruff_cache/` if not covered. Commit before any code change below, and commit each step
separately so the friend's changes land on a clean history.

---

## Step 1 — wire the similarity features into Stage D

All in `src/phageforge/funnel/stages.py`.

**1a. Fix the missing-candidate default.** `_knn_scores` (`stages.py:342`) returns the top
`k` phages per arm; `stage_a_candidates` fills absent ones with `0.0` (`stages.py:502-503`).
ES scores cosine as `(1 + cos) / 2`, so `0.0` means *perfectly anti-similar* — strictly
worse than orthogonal. At 96 phages every candidate is returned so it never bites, but it
is a silent scale bug that would misrank at the tier sizes Step 2 introduces. Default to
`0.5` (orthogonal) instead, and note why.

**1b. Centre the features.** Feed `rbp_similarity - 0.5` and `genome_similarity - 0.5` into
Stage D so a weight of zero is genuinely neutral and a fitted weight is interpretable as
"per unit of cosine above orthogonal". Do the subtraction where the `features` dict is built
(`stages.py:752`), not in the script.

**1c. Populate and consume them.** Add both keys to the per-candidate `features` dict, and
add the two terms to the Painless source (`stages.py:772`). The candidate dicts already
carry the raw values, so this is a lookup, not a new query.

**1d. Keep `FEATURE_NAMES` honest.** It already lists six features (`stages.py:707`); after
1c that tuple finally matches what the script reads. Add a test asserting the two agree —
this is exactly the drift that produced the null ablation.

**Do not hand-pick the two new weights.** Step 1 ships with them still at `0.0` in
`DEFAULT_WEIGHTS`, so behaviour is unchanged and the commit is provably inert. Step 1.5
supplies measured values.

### Step 1.5 — fit the weights on training folds (`bench/weights.py`, new)

`stages.py:692` already promises this: *"`phageforge.bench` fits them on strain-grouped
folds and writes the fitted set to data/derived"*. It is also the step ARCHITECTURE §7.4
wants before the `learning_to_rank` rescorer — same features, same shape, one fewer moving
part.

- Reuse `harness.make_folds`, `harness.load_truth`, `harness.evaluable_strains` and
  `harness.fold_breadth` — do not re-derive the split. **Read-only imports from
  `harness.py`; do not edit it** (peer owns that file).
- For each fold: run the funnel over *training* strains with `top_n=len(phages)` and
  `explain=False, build_cocktail=False, persist=False`, harvest each candidate's `features`
  dict from `run.shortlist`, label from `truth`, fit `sklearn.linear_model.LogisticRegression`
  over `FEATURE_NAMES`. sklearn is already a dependency.
- Write the mean of the per-fold coefficients to `data/derived/weights.json` with the fold
  coefficients alongside, so the variance across folds is visible rather than averaged away.
- Bump `MODEL_VERSION` to `phageforge-linear-v2` when fitted weights become the default.

Then re-run the benchmark and report **both** rows: untrained defaults and fitted weights.
If the RBP arm still contributes nothing once it can actually reach the score, that is a
real, publishable negative result — say so in the README rather than tuning until it moves.

> ⚠️ Re-running `make bench` overwrites `data/derived/benchmark.json`, which peer session
> `0a` wrote at 11:30. Write to `--out data/derived/benchmark-v2.json` first, compare, and
> only promote once the peer's pass has settled.

---

## Step 2 — the hybrid-search benchmark (`bench/hybrid.py`, new)

Answers three questions `current_status.md` raises and never closes.

**Corpus.** There is no real 100K-vector corpus — `pf-proteins` holds 127 docs. Build
`pf-vectorbench` synthetically: resample the 127 real ESM-2 vectors from
`data/derived/cache/esm-esm2_t12_35M_UR50D.npz` and perturb with Gaussian noise calibrated
to the observed pairwise-cosine spread (`cli.py:250` `_cosine_spread` already computes it),
so the corpus has realistic neighbourhood structure rather than uniform noise. **Label every
result "synthetic corpus, ESM-2-derived distribution" in the output and the README.** A
synthetic ANN benchmark is a legitimate measurement of Elasticsearch; it is not a
measurement of phage biology, and the two must not be allowed to blur.

**Tiers: 10K / 50K / 100K.** Not the 500K/1M in the status doc. This box has ~1 GB RAM free
of 7.5 GB with the ES heap already pinned at 1 GB (`infra/elasticsearch.sh`); 500K float32
× 480 dims is 960 MB of vectors and would thrash. Put higher tiers behind `--tiers` and
state plainly in the README that they need the Elastic Cloud move the status doc §8 already
anticipates.

**Measurements**, per tier:

1. **recall@10 of HNSW vs exact.** Ground truth from an exact `script_score` +
   `cosineSimilarity` scan (or numpy over the generated matrix — cheaper and exact).
2. **Latency**, median and p95, over ≥200 queries, after a warm-up pass that is discarded.
3. **Quantization: `hnsw` float32 vs `int8_hnsw` vs `bbq_hnsw`** — recall, p95, and index
   size from `_cat/indices`. This is the §4.1 loop closed: quantization wrecked the
   phylogeny at 404×256, and §4.1 asserts without evidence that it pays at scale. Note that
   `pf-proteins` already ships `bbq_hnsw` (`mappings.py:75`) at 127 docs, far below where
   BBQ is meaningful — if this measurement says so, that is a finding about the current
   mapping, not just about the tiers.
4. **RRF vs pure vector.** Same two-arm query as Stage A (kNN + BM25 over a synthetic
   annotation field), scored against the exact-kNN ground truth, so the fusion earns its
   place or does not.

Reuse `harness.format_table`-style plain-text output and `harness._finite` + `harness.save`
for JSON (import them; do not copy). Write to `data/derived/hybrid.json`. Clean up
`pf-vectorbench` on exit unless `--keep`.

---

## Step 3 — docs

**`README.md`** (new). Short. What it is; the one-paragraph claim and the explicit
non-claims from SOLUTION_OVERVIEW §3; the cold-start command sequence; the benchmark table
as measured; the honest caveats (synthetic hybrid corpus, RBP ablation result whatever it
turns out to be, 390 of 402 strains evaluated because 12 have no positives).

**`current_status.md`** — rewrite §1 (phase table), §2 (the "you cannot run an end-to-end
prediction today" flow diagram is now false), §5.3, and §7. Keep §4's issue log verbatim;
it is the most valuable part of the file. Keep the warning in §3 that the 90% `ECOR-54`
number is not a benchmark — and now point it at the real one.

**CLI + Makefile** (last, re-read before editing): add `phageforge fit-weights` and
`phageforge hybrid` subparsers to `cli.py` mirroring the existing `cmd_bench` shape, and
`fit-weights` / `hybrid` targets to the Makefile.

---

## Verification

```bash
make lint && make test                        # currently 17 ruff errors; peer is on this
make funnel STRAIN=ECOR-54                    # unchanged shortlist after Step 1 (weights still 0)
$(PY) -m phageforge.cli fit-weights           # writes data/derived/weights.json
$(PY) -m phageforge.cli bench --out data/derived/benchmark-v2.json
$(PY) -m phageforge.cli hybrid --tiers 10000  # smoke test before the full run
$(PY) -m phageforge.cli hybrid                # ~10K/50K/100K, writes data/derived/hybrid.json
```

Checks that matter, beyond "it ran":

- **Step 1 is inert.** With weights still `0.0`, `make funnel STRAIN=ECOR-54` must return
  the identical shortlist and scores as before the change. If it moves, 1a/1b are wrong.
- **The ablation now has a channel.** With fitted weights, `funnel` and `funnel_minus_rbp`
  must differ by more than float noise — whichever direction. A second null result after
  Step 1 means the RBP signal is genuinely absent, which is a finding; a null result *before*
  Step 1 meant only that the wiring was broken.
- **No leakage in the fit.** Assert in `tests/test_stages.py` that no strain in a fold's
  held-out set appears in that fold's training strain list.
- **Hybrid recall is sane.** Exact kNN must score recall@10 = 1.0 against itself; if the
  harness cannot reproduce that, the ground truth is wrong and every other number is noise.
- **Watch memory** during the 100K tier (`free -g`, `podman stats phageforge-es`). If the
  cluster goes yellow or the OOM killer fires, drop the top tier and record that the limit
  was hit — the honest result, not a smaller silent one.
