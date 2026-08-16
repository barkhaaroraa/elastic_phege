# PhageForge — Current Status

**Updated:** 2026-08-16 · **Phase:** 5 complete; 6 half-done (quality benchmark ✅, scaling benchmark ⬜)

Working notes on what is built, what actually runs, what broke along the way, and
what is still open. Design intent lives in [`docs/`](./docs/); ground truth about
the dataset lives in [`docs/DATA_NOTES.md`](./docs/DATA_NOTES.md). **Where this
file and the architecture doc disagree, this file is what is true.**

**Picking up work?** Read [`docs/REVIEW_BACKLOG.md`](./docs/REVIEW_BACKLOG.md)
first. It carries the 2026-08-16 domain review, what the benchmark numbers
actually prove, and the workstream split with a file-collision matrix for
parallel sessions.

---

## 1. Where we are

| Phase | State | Evidence |
|---|:--:|---|
| 0 — Foundations | ✅ | ES 9.5.1 green; Python 3.12.14 venv; CPU-only torch |
| 1 — Data | ✅ | Archive MD5-verified, security-audited, parsed |
| 2 — Index layer | ✅ | 5 indices live; ingest idempotent across re-runs |
| 3 — Features | ✅ | Genome vectors validated; RBPs tiered, embedded, denormalised |
| 4 — Funnel | ✅ | N→B→A→D→C orchestrated; runs persisted to `pf-predictions` |
| 5 — API | ✅ | FastAPI over the funnel; 13 tests green; `make lint` clean |
| 6 — Benchmark | 🟡 | **Quality half done and measured.** Scaling half not started |

**An end-to-end prediction now works**, over the CLI and over HTTP.

### What is in Elasticsearch right now

```
pf-interactions   38,435 docs   7,976 infects=true · 30,459 infects=false
pf-bacteria          403 docs   402 with genome_vector
pf-phages             96 docs    96 with genome_vector · 89 with rbp_match_vector
pf-proteins          127 docs   108 high-confidence · 19 heuristic
pf-predictions        12 docs   audit trail, one doc per funnel run
```

All interactions are `evidence_tier: 1`, `resolution: strain`. This is the entire
Tier-1 foundation and the only data anywhere with **real measured negatives**.

---

## 2. Current flow — what actually happens today

```
                    ┌─────────────────────────────────────────┐
   Zenodo 13831957  │  make ingest   (idempotent, ~30 s)      │
   146.8 MB, MD5 ✅ │  matrix + metadata + features → ES      │
                    └────────────────────┬────────────────────┘
                                         ▼
                    ┌─────────────────────────────────────────┐
                    │  make features (~20 s)  make proteins   │
                    │  · phage TNF sketches      → 96/96      │
                    │  · bacteria MDS embeddings → 402/403    │
                    │  · ESM-2 RBP embeddings    → 127 → 89   │
                    └────────────────────┬────────────────────┘
                                         ▼
    ┌────────────────────────── THE FUNNEL ────────────────────────────┐
    │  N  neighbours    ✅  kNN over strain genome vectors             │
    │  B  prior         ✅  similarity-weighted terms agg              │
    │  A  candidates    ✅  RRF over genome + RBP arms — but see §5.3  │
    │  D  ranking       ✅  linear blend via script_score, in-cluster  │
    │  C  explanation   ✅  significant_terms → reason strings         │
    │  F  cocktail      ✅  greedy set-cover on RBP diversity          │
    └──────────────────────────────────────────────────────────────────┘
                                         ▼
                    ┌─────────────────────────────────────────┐
                    │  make api  →  GET /predict/{strain_id}  │
                    │  ~430 ms cold, ~55 ms median in-bench   │
                    └─────────────────────────────────────────┘
```

### Phase 5 — the API

`make api` serves `phageforge.api.main:app`. Interactive docs at `/docs`.

| Endpoint | Purpose |
|---|---|
| `GET /health` | Cluster + per-index counts + **`ready`**: can the funnel actually answer? |
| `GET /predict/{strain_id}` | The funnel. `top_n`, `neighbours`, `explain`, `cocktail`, `persist` |
| `GET /strains` · `/strains/{id}` | Catalogue and type-ahead |
| `GET /phages` · `/phages/{id}` | Catalogue |
| `GET /runs` · `/runs/{run_id}` | The audit trail, summary and full |
| `GET /benchmark?threshold=` | A measured benchmark run. `0` / `1` / `2` select the §3.1 sweep artefacts |
| `GET /` | The demo UI (below) |

### The demo UI

`make api`, then open <http://localhost:8000/>. One self-contained HTML file
served same-origin — no CDN, no build step, no external request (asserted by a
test). It is a thin view over the endpoints above; every number on the page comes
from one of them and nothing is recomputed client-side.

Three panels, matching the three questions someone actually asks:

1. **Pick a strain** — type-ahead over `/strains`, with the strain's receptor
   profile (ST, phylogroup, O-antigen, H-type, LPS core, defence systems) shown
   before you run anything, so the reasons afterwards can be checked against it.
2. **What came back** — the funnel's per-stage latency as a proportional strip
   (which stage actually costs the time), then the ranked shortlist: score,
   evidence tier, neighbour support, and the Stage C reason strings on every row.
   Then the cocktail with its per-phage rationale.
3. **How good is it** — the benchmark against all three baselines and both
   ablations, with the base rate drawn as a reference line, and tabs for the
   three threshold cuts from §3.1. The `− RBP arm` bar sits *above* the funnel
   bar at two of the three cuts; the chart is not arranged to hide that.

`/#ECOR-54` deep-links straight to a result, so a demo link lands on output
rather than an empty form. A strain with no genome vector (`LF110`) renders the
422 stage-N stop as a stated cause rather than an empty list — the designed
behaviour is visible in the UI, not just in the API contract.

Two behaviours are deliberate and tested:

- **An unknown strain is a 404** (with a did-you-mean), **a stopped funnel is a
  422** naming the stage and the cause. A caller can tell "you asked for
  something that does not exist" from "the funnel could not answer", and neither
  is ever a 200 with an empty list.
- **The ablation switches are exposed over HTTP** (`use_rbp`, `use_prior`). They
  are the same code path the benchmark measures, so a reviewer can reproduce an
  ablation row themselves rather than taking the benchmark's word for it.

---

## 3. Benchmark results — the deliverable

5-fold **strain-grouped** CV, 390 evaluable strains × 96 phages, base rate 21.4%,
seed 20260815. A strain appears in exactly one fold, and every other strain in
its fold is excluded from its neighbour set, so nothing leaks.

| Method | P@10 | R@10 | AUPRC | × base | p95 ms |
|---|--:|--:|--:|--:|--:|
| random | 0.204 | 0.098 | 0.246 | 0.95 | — |
| most-generalist | 0.581 | 0.379 | 0.594 | 2.72 | — |
| phylogenetic-NN | 0.585 | 0.385 | 0.590 | 2.74 | — |
| **PhageForge funnel** | **0.633** | **0.413** | **0.654** | **2.96** | 74 |
| funnel − RBP arm | 0.633 | 0.413 | 0.657 | 2.96 | 64 |
| funnel − prior | 0.535 | 0.298 | 0.525 | 2.50 | 65 |

Cocktail of 4 covers **92.1%** of held-out strains, vs 92.3% for the top 4 ranked
individually.

**Readings, including the unflattering one:**

1. The funnel beats both non-trivial baselines — +5.2 pts P@10 over
   most-generalist, +4.8 over phylogenetic-NN. Modest, but it is measured on
   grouped folds and it is real.
2. **The prior is doing the work.** Removing it costs 9.8 pts P@10 and 13 pts
   AUPRC. Collaborative filtering over measured neighbours is the engine.
3. **The RBP arm contributes nothing measurable.** Removing it changes P@10 by
   −0.0003 and *improves* AUPRC by 0.003 — both well inside noise — while saving
   ~6 ms. At 96 phages there is nothing for Stage A to narrow, so this is the
   expected result rather than a surprise, but it must be reported as what it is:
   at this scale the RBP arm is unpaid overhead. §5.3 is where it earns its keep.
4. Most-generalist at P@10 = 0.58 is a reminder of how strong the dumb baseline
   is here. Any claim about this system that is not measured against it is empty.

> The old §3 hand-picked `ECOR-54` figure of P@10 = 90% is **superseded**. It was
> one strain in a dense ST73 cluster. The number to quote is **0.633**.

### 3.1 Threshold sensitivity — the sweep

§5.2's open question, run rather than argued. Each row is a **complete
re-derivation** (`ingest → features → proteins → bench`), not just a relabelled
rerun. Artefacts: `benchmark.json` (`>0`), `benchmark-gt1.json` (`>1` = ≥2),
`benchmark-gt2.json` (`>2` = ≥3).

| Cut | Positives | Base rate | Strains | random | generalist | phylo-NN | **funnel** | − prior | − RBP | **× base** | Cocktail |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| `>0` (default) | 7,976 (20.8%) | 21.4% | 390 | 0.204 | 0.581 | 0.585 | **0.633** | 0.535 | 0.633 | **2.96** | 92.1% |
| `>1` (≥2) | 4,962 (12.9%) | 13.8% | 375 | 0.131 | 0.413 | 0.434 | **0.486** | 0.426 | 0.484 | **3.51** | 85.1% |
| `>2` (≥3) | 2,656 (6.9%) | 8.3% | 335 | 0.078 | 0.274 | 0.298 | **0.341** | 0.272 | 0.347 | **4.11** | 71.0% |

All figures P@10. Fewer strains are evaluable at stricter cuts because a strain
with no positive contributes nothing to P@10 or AUPRC.

**What the sweep settles:**

1. **The ordering never changes.** funnel > phylo-NN > generalist > random at
   every cut. The central claim does not depend on the threshold choice.
2. **Enrichment over base rate improves monotonically** — 2.96 → 3.51 → 4.11.
   Absolute P@10 falls (mechanically, there are fewer positives to find), but the
   funnel gets *relatively better* at surfacing strong lysis. If anything the
   default `>0` **understates** the system.
3. **The prior is the engine at every cut.** Removing it costs 9.8 / 6.0 / 6.9
   pts P@10. At `>2` the ablated funnel (0.272) falls to the generalist baseline
   (0.274) — without neighbour transfer there is nothing strain-specific left.
4. **The RBP arm is null at every cut**, and at `>2` it is *negative* (0.347
   without vs 0.341 with). Three independent thresholds agreeing is much stronger
   evidence than the single run in §3 reading 3.
5. **Cocktail coverage degrades as the bar rises** (92% → 85% → 71%), which is
   the honest reading: guaranteeing *strong* lysis across a strain panel is a
   materially harder problem than guaranteeing *any* lysis, and a 4-phage
   cocktail is not sufficient for it.

**Reproducibility note.** For a fixed index the benchmark is bit-identical across
runs (verified). Across a **re-ingest** it is not: P@10 drifts by up to ~0.2 pts
on the ablation rows, because re-indexing changes segment layout and the
approximate HNSW kNN returns slightly different neighbour sets. Two orders of
magnitude below the effects above, so it changes nothing — but "identical
numbers" should not be promised across a rebuild, only across a rerun.

### Earlier validation (still holds)

| Check | Result | Reading |
|---|---|---|
| TNF nearest-neighbour genus agreement | 91/94 = 96.8% | Composition sketches recover taxonomy |
| MDS reconstruction vs PanACoTA distances | r = 1.0000, rel. MAE 0.08% | Embedding reproduces the phylogeny essentially exactly |
| MDS variance in 256 dims | 100%, 0.04% negative eigenmass | 256 dims loses nothing |
| kNN neighbour sanity (`ECOR-54`) | 8/8 same phylogroup, 7/8 same ST | Nearest 3 share ST73 **and** O25 |
| Stage C on `T4LD` | `lps_type: R1` top signal, 69% vs 51% background | **Recovered known T4 receptor biology unprompted** |

---

## 4. Issues hit, and how they were resolved

### 4.1 int8 quantization silently corrupted the phylogeny — **fixed**

A strain failed to retrieve **itself** as its own nearest neighbour: implied
self-distance 0.0043 instead of 0, ranked 6th.

MDS coordinates decay sharply in magnitude across dimensions; int8 uses a single
scale per vector, so trailing dimensions quantized to zero and the fine structure
was destroyed. Quantization error (~4e-3) exceeded typical inter-strain distance.

**Fix:** plain `hnsw`, float32. At 404 × 256 that costs 413 KB — the saving was
never worth it. *Quantization is the right call at the 100K–1M benchmark tiers,
which makes it a variable to **measure** there rather than assume.*

### 4.2 `l2_norm` score resolution exhausted — **fixed**

PanACoTA distances are tiny (max 0.117, typically ~3e-4). Elasticsearch scores
`l2_norm` as `1/(1+d²)`, so every score landed within 1e-7 of 1.0 — beyond
float32's ability to distinguish, which would have made Stage B's similarity
weighting meaningless.

**Fix:** scale MDS coordinates to unit RMS (×29.4). A uniform similarity
transform — ranking is identical, numeric range is restored.

### 4.3 ~2 GB of CUDA libraries on a GPU-less machine — **fixed**

The default PyPI `torch` wheel pulled `nvidia-cusolver`, `nvidia-nccl`,
`nvidia-cuda-cupti` and friends. 1.6 GB was cached before it was caught.

**Fix:** pinned torch to the CPU index in `pyproject.toml`. 502 MB + ~2 GB → **183 MB**.

### 4.4 Background downloads silently killed — **fixed (process hygiene)**

Wrapping `nohup … &` *inside* a backgrounded tool call meant the wrapper exited
immediately, the harness reported success, and the process group was reaped —
killing the real download mid-transfer. Cost two false starts.

**Fix:** never nest the two. Background the command directly.

### 4.5 ES 9.x omits `dense_vector` from `_source` — **behaviour, not a bug**

Vectors are excluded from the default `_source` response as a bandwidth
optimisation and returned only when explicitly requested via `source_includes`.
Naive code sees no vector and fails confusingly. Worth remembering — the API's
`has_rbp_vector` flag infers presence from `rbp_match_id` for exactly this reason.

### 4.6 `benchmark.json` was not valid JSON — **fixed**

Metrics that genuinely do not apply to a method (cocktail coverage for `random`,
latency for offline baselines) are carried as `nan`, which is correct internally.
But `json.dumps` writes that as a bare `NaN` token — valid Python, **invalid
JSON**. `jq`, `JSON.parse` and any strict parser rejected the entire file, which
would have broken the `/benchmark` endpoint and any downstream chart.

**Fix:** `harness._finite()` converts non-finite floats to `null` before writing,
and `save()` now passes `allow_nan=False` so a regression fails loudly at write
time rather than silently producing an unreadable artefact. The existing results
file was rewritten in place; the measured numbers are unchanged.

---

## 5. Open issues

### 5.1 ✅ RBP ordinal mapping — **decided and shipped**

Resolved as proposed: confidence-tier rather than guess, and never average the
tiers away. Of 130 identifiers in `RBP_list.csv`:

| Tier | Count | Treatment |
|---|--:|---|
| `high` — exact ordinal, plausible length | 108 | embedded |
| `heuristic` — recovered within ±3, flagged | 19 | embedded, flagged |
| `unresolved` — ordinal past end of genome | 3 | skipped |

127 proteins embedded with ESM-2 `esm2_t12_35M_UR50D` (480-dim, CPU), centred
against the corpus mean, and denormalised onto **89 of 96 phages** as
`rbp_match_vector`. The tier travels with the data: `pf-proteins.confidence`,
`pf-phages.rbp_match_confidence`, and out through the API on every shortlist
entry and cocktail member. The ±3 window was *not* applied blindly — each
heuristic match carries a `resolution_note` a reviewer can read.

Running Prokka to reproduce the authors' numbering exactly remains the clean
fix. Still deferred, still a heavy dependency, and now demonstrably low-value
(§3, reading 3).

### 5.2 🟢 Infection threshold — **swept and measured; the reviewer's call is now low-risk**

The matrix is a graded 0–4 lysis score, not binary and not EOP:

| Score | 0 | 1 | 2 | 3 | 4 | empty |
|---|--:|--:|--:|--:|--:|--:|
| Cells | 30,459 | 3,014 | 2,306 | 1,368 | 1,288 | 157 |

Default is `> 0` (7,976 positives, 20.8%), grounded in the authors' own code,
which compares `>0` throughout. **This is review question 3 in
`SOLUTION_OVERVIEW.md` §12.** It was flagged as the largest uncertainty under
every number in §3. It has now been run at all three defensible cuts — see
[§3.1](#31-threshold-sensitivity--the-sweep). **The conclusions are
threshold-invariant**: the method ordering is identical at every cut and the
funnel's enrichment over base rate *improves* as the cut tightens.

The reviewer still chooses which cut to publish. What has changed is that the
choice no longer threatens the headline claim — it rescales the absolute
numbers, and the ranking survives.

> ⚠️ **Off-by-one trap.** The code applies `score > PF_INFECT_THRESHOLD`, so the
> "≥2" cut this section asks about is `PF_INFECT_THRESHOLD=1`, **not** `=2`
> (which is ≥3). Setting `=2` and labelling it "≥2" silently reports the wrong
> experiment.

> ⚠️ **A threshold change is not a one-command rerun.** `infects` is computed at
> *ingest*, and `bulk_index` replaces whole documents — so a re-ingest wipes the
> `genome_vector` and `rbp_match_vector` fields that `features` and `proteins`
> attach afterwards. The full chain is
> `ingest → features → proteins → bench` (~3 min; ESM-2 is served from cache).
> Running `bench` alone with a different env var produces a file *labelled* with
> the new threshold while measuring the old data.

### 5.3 🟠 Stage A still has nothing to do at this scale — **now measured**

No longer a suspicion: §3 reading 3 shows the RBP arm changes nothing at 96
phages. Narrowing 96 candidates to 500 is a no-op by construction. Stage A only
becomes real at the hybrid-search scale tiers, which is precisely the half of
Phase 6 that has not been built (§5.6).

### 5.4 🟡 Small data-coverage gaps

- `LF110` is indexed but absent from the distance matrix → no `genome_vector`
  (402/403). `/health` reports it; `/predict/LF110` correctly 422s naming stage N.
- `H1-005-0065-L-P` and `H27` are in the matrix but not in `pf-bacteria`.
- 1 of 97 phage FASTAs has no matching indexed phage.
- 7 of 96 phages have no usable RBP and so no `rbp_match_vector`; the cocktail
  labels these "receptor diversity unverified" rather than assuming no overlap.
- Our counts are **38,435 measured**, not the paper's 38,688 (403 × 96). One
  strain has metadata but no interaction row. Report ours.
- The benchmark evaluates **390 of 402** strains; the remainder have too few
  measured phages to score. Reported in the output, not hidden.

### 5.5 🟡 Deferred, with re-entry points

Bacterial genomes and protein sequences are **not in the archive**, so:
`lifestyle` is `unknown` for all phages (BACPHLIP not run); there are no
bacterial receptor vectors, so Stage A's third arm is RBP-similarity to phages
infecting nearby strains rather than a direct biophysical match; cocktail design
uses RBP-embedding cosine as its diversity proxy rather than derived
`receptor_classes`.

### 5.6 🔴 Phase 6's scaling half is not built — **the main gap**

The quality benchmark is done. The hybrid-search benchmark is not: recall@10 vs
exact kNN, p95 latency, and RRF-vs-pure-vector across 100K / 500K / 1M vector
tiers. That is where quantization (§4.1), Stage A (§5.3) and the move to Elastic
Cloud all become measurable rather than assumed — and on a 7.5 GB box with a 1 GB
ES heap, the 500K and 1M tiers will not run locally.

### 5.7 ⚪ Non-technical

Hackathon deadline unconfirmed — public listings say submissions closed
**13 Aug 2026**, with Elastic{ON} Mumbai 30 Sept for shortlisted teams. Worth
checking the registration email; it changes the milestone maths, not the code.

---

## 6. Architecture deviations from `docs/ARCHITECTURE.md`

Five, all forced by what the data actually contains:

1. **Bacteria `genome_vector` is MDS over phylogenetic distances**, not k-mer
   sketches — no bacterial genomes are shipped. Arguably better: kNN reproduces
   true phylogenetic neighbourhood rather than approximating it.
2. **Bacteria vectors use `l2_norm`; phage vectors use `cosine`.** MDS
   coordinates are Euclidean by construction; TNF profiles are compositional.
   Using cosine on MDS coordinates would compare directions from an arbitrary
   origin and distort the phylogeny.
3. **No quantization at this scale** (§4.1).
4. **Stage A's RBP arm is reshaped** (§5.5) and currently inert (§5.3).
5. **Stage D is a linear blend in `script_score`, not an LTR rescorer.** Same
   shape one step earlier: features come from the previous stages and the
   *scoring runs inside Elasticsearch*. Swapping in the XGBRanker via Eland is a
   change of one query body, not of the architecture — and it needs the trained
   model that Phase 6 was going to produce.

---

## 7. Next

1. **Re-run the benchmark at `PF_INFECT_THRESHOLD=2`** (§5.2) and report both.
   One command; it is the largest single uncertainty in every number above.
2. **Build Phase 6's scaling half** (§5.6) — the tiers, and with them the case
   for Stage A and for quantization. Needs Elastic Cloud beyond 100K.
3. **Train the LTR model** on the fold structure the harness already produces and
   swap Stage D's `script_score` for the rescorer.
4. Optional: BACPHLIP for `lifestyle`, if bacterial genomes can be sourced.

**Phase 6's scaling half is what is left of the deliverable.**

---

## 8. Reproducing from scratch

```bash
make venv          # uv + Python 3.12 venv, CPU-only torch
make up            # Elasticsearch 9.5.1 via podman, waits for green
make ingest        # verify MD5 → extract → parse → index (idempotent)
make features      # TNF sketches + MDS embeddings
make proteins      # RBP recovery (pyrodigal) + ESM-2 embeddings
make info          # document counts per index
make bench         # the benchmark (~70 s)
make api           # FastAPI on :8000, docs at /docs
make test lint     # 13 tests, ruff clean
```

Environment: Python 3.14 is the only system Python and is too new for the bio
stack, hence the pinned 3.12 venv. No Docker — **podman**. No GPU. 7.5 GB RAM
total, which is why the ES heap is capped at 1 GB and why the move to Elastic
Cloud is expected at Phase 6's scale tiers.
