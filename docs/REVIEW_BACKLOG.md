# PhageForge — Review Findings & Work Backlog

**Created:** 2026-08-16 · **Purpose:** hand-off document for parallel sessions

This file exists so that any session — or any person — can pick up work without
re-deriving the state of the project. It records (1) what is *verifiably* true
right now, (2) what a domain reviewer found, (3) what the benchmark actually
proves, and (4) discrete workstreams scoped so that two sessions working at once
do not collide.

> **Precedence.** `current_status.md` (root) is the authority on **what is
> built**; it was refreshed 2026-08-16 11:36 and is current. This file is its
> complement: **what the numbers mean, what the reviewer found, and who does
> what next.** It does not restate build state except where a finding depends on
> it. `docs/ARCHITECTURE.md` is design *intent*, not state — where it disagrees
> with either of the other two, it is the one that is out of date.
> Verify before trusting any of them (§1.1).

---

## 1. Verified state — 2026-08-16

### 1.1 How to re-verify in 30 seconds

Do not trust this section. Re-run it — the project moves fast.

```bash
make status                                   # cluster health
curl -s 'localhost:9200/_cat/indices/pf-*?v'  # document counts
.venv/bin/pytest -q && .venv/bin/ruff check src tests
python -c "import json;print(json.load(open('data/derived/benchmark.json'))['results'])"
```

### 1.2 Indices

| Index | Docs | Note |
|---|--:|---|
| `pf-interactions` | 38,435 | 7,976 positive · 30,459 negative · all Tier 1, strain-resolution |
| `pf-bacteria` | 403 | 402 with `genome_vector` (`LF110` has none) |
| `pf-phages` | 96 | all with `genome_vector` |
| `pf-proteins` | 127 | ESM-2 RBP embeddings; 108 high-confidence · 19 heuristic |
| `pf-phages` (RBP) | 89 | denormalised `rbp_match_vector` — 7 phages have none |
| `pf-predictions` | 12 | audit trail, written per funnel run |
| `pf-literature` | — | **does not exist**; needs ELSER on Elastic Cloud |

### 1.3 Code

All of Phases 3–6 exist and run. `pytest` = **13 passed**, `ruff` = clean.

| Module | Lines | State |
|---|--:|---|
| `funnel/stages.py` | 847 | Stages A/B/C/D all implemented |
| `funnel/pipeline.py` | 257 | orchestration, `FunnelRun`, persistence |
| `funnel/cocktail.py` | 137 | greedy set-cover |
| `api/main.py` | 454 | 9 endpoints, incl. `GET /predict/{strain_id}` with ablation switches |
| `bench/harness.py` | 551 | strain-grouped folds, 3 baselines, 2 ablations |
| `features/proteins.py` | 517 | RBP recovery + ESM-2 |

Test coverage is **API-only** (`tests/test_api.py` is the sole test file).
Nothing covers the funnel stages, features, ingest, or the harness.

### 1.4 Benchmark — `data/derived/benchmark.json`

390 strains · 5 folds · 96 phages · `infect_threshold = 0.0` · base rate 0.2139
· seed 20260815 · 69.3 s wall.

| Method | P@10 | R@10 | AUPRC | p95 ms |
|---|--:|--:|--:|--:|
| random | 0.204 | 0.098 | 0.246 | — |
| generalist | 0.581 | 0.379 | 0.594 | — |
| phylo_nn | 0.585 | 0.385 | 0.590 | — |
| **funnel** | **0.633** | **0.413** | **0.654** | 73.7 |
| funnel − rbp | 0.633 | 0.413 | 0.657 | 64.3 |
| funnel − prior | 0.535 | 0.298 | 0.525 | 65.2 |

Funnel cocktail coverage 0.921 · top-4 coverage 0.923 · 0 failures.

---

## 2. What the benchmark actually proves — read this before claiming anything

Three readings, all uncomfortable, all supported by the table above.

### 2.1 🔴 We are +4.8 points over a phylogenetic nearest-neighbour lookup

`phylo_nn` — "find the nearest strains by phylogeny, recommend their phages" —
scores 0.585. The full funnel scores 0.633. Every stage, vector, and aggregation
buys **4.8 points of P@10** over a baseline that fits in twenty lines.

Against `generalist` ("rank phages by how many strains they infect, ignore the
query entirely") the margin is +5.2 points.

This is the number a reviewer or judge will press on. Until §3.1 lands,
"strain-level and explainable" is a **design claim, not a demonstrated
advantage**. Do not present the funnel as validated without stating this margin.

### 2.2 🔴 The RBP arm contributes exactly nothing

`funnel_minus_rbp` scores **0.6333** vs the full funnel's **0.6331** — the
ablation is fractionally *better*, and its AUPRC is higher too (0.657 vs 0.654).
It is also 9 ms faster at p95.

All of Phase 3 — the risky ordinal mapping (§5.1 of `current_status.md`), the
127 ESM-2 embeddings, the confidence tiering — currently buys zero. §5.3
predicted this ("Stage A has nothing to do at this scale"); it is now measured.

### 2.3 🟠 The prior is the only thing doing real work

`funnel_minus_prior` collapses to 0.535 — below both non-random baselines. Stage
B's neighbour-transfer prior *is* the system. Everything else is decoration on
top of it, which is consistent with 2.1 and 2.2.

---

## 3. Domain review — 2026-08-16

Reviewed by **Barkha Arora** (microbiology). This partially satisfies the
"reviewed by a domain expert" line in `TECHNICAL_DESIGN.md` §11 and answers
several of the questions posed in `SOLUTION_OVERVIEW.md` §12.

### 3.1 🔴 CONFIRMED — genetic similarity does not imply capsule similarity

> *"genetic similarity doesn't mean capsules of bacteria will also be similar;
> sometimes it is the deciding factor for infections and resistance … it may just
> pick up bacterias/phages that are almost exactly copies in terms of dna but
> that doesn't mean it will function the same."*

**This is correct, and it is a structural defect in Stage D, not a data gap.**

We *have* the receptor data and we *index* it:

| Field | Values present |
|---|---|
| `lps_type` | R1 (206) · R3 (61) · K12 (43) · No_waaL (35) · R4 (35) · R2 (22) |
| `capsule_types` | GroupIV_s (353) · Wzy_stricte (312) · GroupIV_e_stricte (218) · ABC (142) |
| `o_antigen` | Oneg (40) · O2 (29) · O8 (25) · O6 (21) · O25 (19) · O9 (18) … |

But every one of these appears **only** in `EXPLAIN_FIELDS`
(`src/phageforge/funnel/stages.py:526`) — they are consumed by Stage C to
generate reason strings, and by nothing else.

The Stage D scoring features (`stages.py:694`) are:

```python
DEFAULT_WEIGHTS = {
    "prior": 3.0, "support": 0.30, "candidate_rank": 0.50, "breadth": 0.25,
    "rbp_similarity": 0.0,      # declared, weighted zero
    "genome_similarity": 0.0,   # declared, weighted zero
    "bias": 0.0,
}
```

**So capsule and LPS currently write the sentence explaining the ranking without
influencing the ranking.** That is the reviewer's critique, exactly, located to a
line number. It also explains §2.1: with receptor biology excluded from scoring,
there is little left to beat phylogeny with.

Corroborating evidence that the signal is real and being wasted: Stage C already
surfaced `lps_type: R1` as the top signal for phage T4LD at 69% vs 51%
background — recovering known T4 receptor biology unprompted.

### 3.2 🟠 CONFIRMED with a correction — CRISPR / existing resistance

> *"some bacterial strains will already have resistance to a phage … CRISPR
> sequences … you'll have to incorporate that so it checks existing resistances."*

The premise is right and we are **better positioned than the reviewer assumed** —
no external CRISPR-finding tool is needed. DefenseFinder was run by the original
authors across all 403 strains and the output is already indexed in
`pf-bacteria.defence_systems`:

| System | Strains |
|---|--:|
| `MazEF` | 338 |
| `RM_Type_I` | 327 |
| `RM_Type_IV` | 268 |
| **`CAS_Class1-Subtype-I-E`** | **215** |
| `Mok_Hok_Sok` | 193 |
| `RM_Type_II` | 155 |
| `CAS_Class1-Subtype-I-F` | 66 |

…plus ~100 more (BREX, CBASS, Gabija, Zorya, Thoeris, Septu, Druantia…).

**Two corrections worth carrying into any write-up:**

1. **CRISPR presence ≠ resistance to a given phage.** Resistance requires a
   *spacer matching that specific phage's genome*. We have system-level
   detection, not spacer sequences. Real spacer matching needs the bacterial
   genomes, which are **not in the Zenodo archive** — this is the concrete
   reason to pull genomes from NCBI/Enterobase (WS-2).
2. **In this dataset restriction–modification systems dominate CRISPR** (327 vs
   215). The anti-phage defence story here is broader than CRISPR alone, and a
   `defence_systems` scoring feature should cover all of it, not just Cas.

Like the receptor fields, `defence_systems` is currently explanation-only.

### 3.3 🟡 PARTIALLY ACCEPTED — docking / UniProt / PDB / AutoDock Vina

> *"that will need docking simulations … you need to add in uniprot/pdb data …
> autodock vina is what's mainly used … I think u should add the docking/protein
> interaction part into the scoring system."*

The *direction* is right and is reinforced by §2.2 — the sequence-level protein
arm is inert, so structure is the obvious next lever. But the scale must be
stated honestly rather than promised:

- **Not feasible:** pairwise docking across 96 phages × 403 strains ≈ 38,000
  Vina runs, each needing structures on both sides. This is out of reach on the
  current hardware (no GPU, 7.5 GB RAM) and out of scope for the timeline.
- **Feasible:** ESMFold over the **127 RBPs we already have** — small model, CPU,
  hours not weeks.
- **Feasible as a rescorer:** Vina on the **top ~10 candidates per strain only**,
  as a secondary confidence indicator layered onto an existing ranking — never as
  a primary filter.

Sequencing note: **§3.1 must be tried first.** It is cheaper by an order of
magnitude, uses data already in hand, and tests the same underlying hypothesis
(that receptor compatibility beats phylogeny). If receptor *features* do not
move the number, receptor *structures* are a much larger bet on the same idea.

### 3.4 🟡 ACCEPTED — dataset is E. coli only

> *"the dataset you're using is only for ecoli and its phages … there's a similar
> pipeline for klebsiella already but that's only for klebsiella."*

Correct, and already planned: PhageHostLearn / Zenodo 11061100 (WS-7). Its
importance has increased — a second organism is the test that the pipeline
**generalises** rather than fitting Gaborieau's particular matrix.

### 3.5 🟢 ACCEPTED — positioning and prior art

> *"needs to be marketed as prediction + database for pre-wet-lab analyses"*

This is the right frame and largely describes what already exists: evidence
tiers, the `pf-predictions` audit trail, cocktail set-cover, and the planned
96-well plate export. It also defuses the reviewer's adoption concern (phage
therapy being far more established in Georgia/Russia than in India) — as a
research triage tool, it does not depend on clinical phage therapy being
widespread locally.

**Prior art to acknowledge explicitly in any write-up** (the reviewer is right
that we are building on an ecosystem, not founding a field):

| Tool | What it does | How we differ |
|---|---|---|
| HostPhinder (2016) | phage → host *genus/species* via shared k-mers; 81% genus | one phage → many hosts; we do one strain → ranked phages |
| WIsH · Phirbo · VirHostMatcher-Net | k-mer / Markov host prediction | same species-level framing |
| RaFAH · PHERI · vHULK | ML on genomic features | no measured negatives, no strain resolution |
| **PhageHostLearn (2024)** | strain-level RBP↔receptor ML, *Klebsiella*, AUC 81.8%, wet-lab validated | closest prior art; different organism. Gaborieau 2024 is the E. coli analogue we build on |

A 2024 benchmark reviewed **27** such tools and found inconsistent evaluation
across them. Our defensible niche is (a) strain-level rather than species-level,
(b) measured negatives, and (c) retrieval + evidence + inventory in one ranked,
explained output — **not** a novel prediction algorithm. Note also that our
`phylo_nn` baseline is essentially the similarity-transfer logic those tools use,
so we already report ourselves against the family honestly.

---

## 4. Workstreams

Each workstream is independently startable. **Check the collision matrix (§4.9)
before starting** — two sessions editing `stages.py` at once will conflict.

Priority: 🔴 blocks the central claim · 🟠 significant · 🟡 valuable · ⚪ optional.

### WS-1 🔴 Receptor & defence features in ranking

**The single highest-value item in this document.** Directly addresses §3.1 and
§3.2, uses data already indexed, and is the cleanest test of §2.1.

Promote the `EXPLAIN_FIELDS` from explanation to scoring. Two new Stage D
features:

- **`receptor_compat`** — fraction of this phage's known host strains sharing the
  target strain's `lps_type` / `capsule_types` / `o_antigen`. This converts Stage
  C's `significant_terms` output from prose into signal.
- **`defence_mismatch`** — whether the target carries defence systems that are
  *depleted* among this phage's known hosts (covering RM, BREX, CBASS… not only
  Cas, per §3.2).

Both must be computed **fold-safe**: like the existing `breadth` override in
`stage_d_rank`, they must be derived from training strains only, or the benchmark
leaks. Read the `breadth` docstring at `stages.py:719` first — the leak it guards
against is the same one.

- **Files:** `funnel/stages.py` (`DEFAULT_WEIGHTS`, `FEATURE_NAMES`,
  `stage_d_rank`), `bench/harness.py` (fold-safe feature computation)
- **Done when:** benchmark rerun with a new `MODEL_VERSION`, and the result —
  better *or worse* — is recorded in §1.4 and §2.
- **Note:** a null result is a real finding. If receptor features do not beat
  0.633, the reviewer's hypothesis is falsified *in our data*, which is worth
  reporting rather than burying.

### WS-2 🟠 CRISPR spacer resistance

Blocked on data: bacterial genomes are not in the archive. Fetch from
NCBI/Enterobase (anonymous; a free NCBI key raises the rate limit from 3 to 10
req/s), run CRISPRDetect or CRISPRCasFinder, extract spacers, match against the
96 phage genomes we already have. A spacer hit is a strong *negative* signal.

- **Files:** new `ingest/genomes.py`, new `features/crispr.py`, `es/mappings.py`
  (add `crispr_spacers`), then one feature in `stages.py`
- **Depends on:** coordinate with WS-1 before touching `stages.py`

### WS-3 🟠 Hybrid-search benchmark tiers

The half of Phase 6 that does not exist. `benchmark.json` has no retrieval
metrics at all: recall@10 vs exact kNN, p95 latency by corpus size, RRF vs pure
vector, quantization at 100K/500K/1M vectors. This is where Stage A and Elastic
Cloud earn their place — and per §2.2, Stage A currently has no measured value at
n=96, so this is its only route to justification.

Needs the INPHARED corpus (~30k phage genomes, AGPL-3.0, anonymous download) to
reach the scale tiers. Revisit the int8 finding from `current_status.md` §4.1 —
it explicitly deferred quantization to "a variable to measure" at these tiers.

- **Files:** new `bench/hybrid.py`, new `ingest/inphared.py`
- **Collides with:** nothing in WS-1

### WS-4 🟠 LTR rescorer

`stage_d_rank` is an honest linear `script_score` placeholder, deliberately
LTR-shaped — its docstring says swapping in the real rescorer is "a change of one
query body, not of the architecture." Train XGBRanker on strain-grouped folds,
export via Eland to `phageforge-ltr-v1`.

- **Files:** new `bench/trainer.py`, `funnel/stages.py` (`stage_d_rank`)
- **⚠️ Sequence after WS-1** — both rewrite Stage D's feature vector. Training a
  ranker on the current four features before adding the receptor features would
  mean retraining immediately.
- **Risk R5:** `learning_to_rank` is technical preview. Pin the version.

### WS-5 🟡 Literature / ELSER — Stage E

`pf-literature` does not exist and cannot on the local cluster: the
`.elser-2-elasticsearch` inference endpoint is Cloud-only, which is why
`config.LOCAL_INDICES` deliberately excludes it. Requires the Elastic Cloud API
key (`PF_ES_API_KEY`, already wired at `config.py:24`). Source is Europe PMC OA
REST — public, no key, rate-limited.

- **Files:** `es/mappings.py`, new `ingest/literature.py`, `funnel/stages.py`
  (new `stage_e_evidence`), `api/main.py`
- **Blocked on:** Elastic Cloud provisioning

### WS-6 🟡 Serving completeness

Three gaps against `TECHNICAL_DESIGN.md` §11:

- **FASTA upload → prediction.** `/predict` accepts a `strain_id` only; there is
  no POST and no upload path. This is an explicit definition-of-done line.
- **96-well plate CSV** (`ARCHITECTURE.md` §8.4). Not implemented anywhere.
- **Frontend.** Nothing exists. §8.5 calls the funnel visualisation "the demo" —
  10,000 → 500 → 50 → 10 narrowing, each stage clickable.

- **Files:** `api/main.py`, `api/models.py`, new `web/`
- **Collides with:** WS-5 (both touch `api/main.py`)

### WS-7 🟡 Second organism — Klebsiella

Ingest PhageHostLearn / Zenodo 11061100 per §3.4. Reuse their extraction
notebooks; do not reimplement.

- **Files:** new `ingest/phagehostlearn.py`
- **Collides with:** nothing

### WS-8 🟡 Test coverage

Coverage is API-only. Nothing tests the funnel stages, features, ingest, or the
harness — the ablation logic in particular is untested and is now load-bearing
for the project's central claim (§2.2 rests entirely on it).

Start with `bench/harness.py`: fold construction, the `breadth` leak guard, and
the ablation switches. If the ablation logic is wrong, §2 is wrong.

- **Files:** `tests/` only
- **Collides with:** nothing — good first task for an extra session

### WS-9 ⚪ Structural docking

Per §3.3. **Do not start before WS-1 reports.** Scope it as ESMFold over the 127
existing RBPs, then Vina as a top-10 rescorer only.

### 4.9 Collision matrix

| File | WS-1 | WS-2 | WS-3 | WS-4 | WS-5 | WS-6 | WS-7 | WS-8 |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| `funnel/stages.py` | ●● | ● | | ●● | ● | | | |
| `bench/harness.py` | ●● | | ● | ● | | | | |
| `es/mappings.py` | | ● | | | ● | | ● | |
| `api/main.py` | | | | | ● | ●● | | |
| `ingest/*` (new files) | | ● | ● | | ● | | ● | |
| `tests/` | | | | | | | | ●● |

●● = primary owner, do not edit concurrently · ● = touches

**Safe to run fully in parallel:** WS-1 + WS-3 + WS-7 + WS-8.

---

## 5. Open decisions — not ours to make alone

### 5.1 🟢 Infection threshold — **swept 2026-08-16; no longer a project risk**

`benchmark.json` confirms `infect_threshold: 0.0`. The matrix is a graded 0–4
lysis score; binarising at `> 0` gives 7,976 positives (20.8%), grounded in the
original authors' own code. The worry was that if the cut belongs at ≥ 2,
positives fall to 4,962 and **every number in §1.4 moves.**

**The numbers move; the conclusions do not.** All three defensible cuts have now
been run end to end (`current_status.md` §3.1):

| Cut | Positives | Base | funnel P@10 | best baseline | × base |
|---|--:|--:|--:|--:|--:|
| `>0` | 7,976 | 21.4% | 0.633 | 0.585 | 2.96 |
| `>1` (≥2) | 4,962 | 13.8% | 0.486 | 0.434 | 3.51 |
| `>2` (≥3) | 2,656 | 8.3% | 0.341 | 0.298 | 4.11 |

Method ordering is identical at every cut, and enrichment over base rate
*improves* as the cut tightens — the default `>0` understates the system rather
than flattering it. The reviewer still picks the cut to publish; that choice is
now a presentation decision, not a threat to the central claim.

Two corrections to the framing above, both of which would have produced a wrong
experiment:

- The code applies `score > PF_INFECT_THRESHOLD`, so **"≥2" is
  `PF_INFECT_THRESHOLD=1`**, not `=2`.
- **It is not a 69 s rerun.** `infects` is set at ingest and `bulk_index`
  replaces whole documents, so a re-ingest drops the vectors attached later by
  `features` / `proteins`. The chain is `ingest → features → proteins → bench`
  (~3 min). Running `bench` alone with the env var set mislabels the old data.

**Consequence for §6:** "Settle §5.1 with the reviewer, rerun, record" is done
apart from the reviewer's preference, and **WS-1 is no longer gated on it** — a
Stage D rework can proceed against `>0` knowing the ordering holds elsewhere.

### 5.2 🟠 Remaining §12 review questions

§3 above answers the receptor-vocabulary and resistance questions. Still open:
the "hit" definition (5.1), whether the receptor-class vocabulary is right
(risk R7), and the blunt one in §12.

### 5.3 ⚪ Deadline

`current_status.md` §5.6: public listings said submissions closed **13 Aug 2026**
— three days before this file was written. Unconfirmed against the registration
email. It changes milestone maths, not code, and the system stands as a portfolio
artifact regardless.

---

## 6. Suggested order

1. **WS-1** (hours) — the one experiment that could change the project's central
   claim. Everything in §2 argues for doing this before anything else is built on
   top of Stage D.
2. **Settle §5.1** with the reviewer, rerun, record.
3. **WS-3** — the remaining route to justifying Stage A and the Elastic Cloud
   move.
4. Then WS-4 / WS-7 / WS-6 in parallel as capacity allows.

WS-8 is unblocked at any time and is the right pick-up for a spare session.

**Do not build WS-4, WS-5, WS-6 or WS-9 on top of the current Stage D** — WS-1
changes its feature vector, and anything trained or presented against the present
0.633 will need redoing.

## 7. Second-pass review — corrections to this document

**Added 2026-08-16 by a parallel session.** Appended rather than edited in place
because §1–§6 were being written concurrently. Each item names the section it
amends; fold them in when the file is next quiet.

Everything in §1 was re-verified independently and is exact: index counts, the
`lps_type` / `capsule_types` / `o_antigen` / `defence_systems` distributions,
`ruff` clean, 13 tests passing. §5.1's threshold sweep supersedes the open
question it replaced and is correct as written. The corrections below are to §2,
§4 and §3.5.

### 7.1 🔴 Amends §2.2 — the RBP arm is *disconnected*, not *worthless*

§2.2 concludes "the RBP arm contributes exactly nothing… all of Phase 3 currently
buys zero". §3.1 shows why, but does not connect it back:

- `stage_a_candidates` computes per-arm similarity into every candidate dict
  (`stages.py:502-503`).
- `stage_d_rank` never reads it. `FEATURE_NAMES` (`stages.py:707`) declares six
  features; the `features` dict (`stages.py:752`) builds four; the Painless source
  (`stages.py:772`) sums four. `rbp_similarity` and `genome_similarity` are
  weighted `0.0` and referenced nowhere.

So `use_rbp=False` reaches the score only through Stage A's RRF reordering
feeding `candidate_rank`. **The ablation measures the reordering channel, not the
RBP signal.** "Phase 3 buys zero" is a claim about wiring, not about ESM-2.

Restate §2.2 as: *the RBP arm is disconnected from scoring, so its value is
currently unmeasured.* The correct reading of §2.3 is unchanged — the prior is
still the only feature carrying weight — but it is doing so partly by default.

### 7.2 🟠 Amends §1.4 and §2.2 — those numbers are from different index states

§2.2's `0.6333 vs 0.6331` comes from `benchmark-t0.json` (11:44). The current
`benchmark.json` (11:50) reads **0.6341 vs 0.6331** — a 0.0010 gap, 5× the quoted
one. `cmd_proteins` was re-run at 11:49, rewriting `esm-mean.npy` and every
`rbp_match_vector`, which moves all three funnel rows.

Back-to-back runs *are* deterministic — three consecutive `bench` runs restricted
to `phylo_nn,generalist` were bit-identical — so this is not sampling noise. It is
state drift: the archived benchmark files were produced against different index
contents and are not comparable at four decimals.

**Do:** record the index state alongside each benchmark file (a
`_cat/indices` snapshot, or at minimum `pf-proteins` and `pf-phages` doc counts
plus the `esm-mean.npy` mtime), and quote no more than three decimals across
files.

### 7.3 🔴 New — there are no error bars anywhere, and per-strain scores are discarded

`MethodScores` holds `p_at_k` as a per-strain list (390 paired observations) and
`summary()` collapses it to a mean. Only the mean is saved, so **significance
cannot be computed after the fact.**

Without a paired test, both "+4.8 points" (§2.1) and "exactly nothing" (§2.2) are
assertions. A Wilcoxon signed-rank or paired bootstrap over the lists already in
memory is ~15 lines.

This matters most for WS-1: a receptor feature worth 0.005 P@10 is real and
useful, and indistinguishable from nothing without a confidence interval.
**Persisting per-strain scores and adding a paired test belongs in §6 ahead of
WS-1**, not in WS-8.

### 7.4 🔴 Amends WS-1 — it will reproduce the §2.2 bug unless plumbing is fixed first

WS-1 adds `receptor_compat` and `defence_mismatch` to `DEFAULT_WEIGHTS` and
`FEATURE_NAMES`. That is precisely the pattern that left `rbp_similarity` dead.
Added the same way without touching the Painless source, they become two more
silently inert features and another meaningless null result.

**WS-1 step zero:** make the script consume `FEATURE_NAMES`, and add a test
asserting the two agree. Also fix `stages.py:502` — absent candidates default to
`0.0`, but ES scores cosine as `(1 + cos) / 2`, so `0.0` means *anti-similar*,
not *unknown*; the neutral value is `0.5`. Harmless at n=96 where Stage A returns
the whole corpus, a silent misranking at WS-3's tier sizes.

### 7.5 🟠 Amends §4.9 — the matrix contradicts its own conclusion

The matrix marks `bench/harness.py` as WS-1 primary (●●) *and* WS-3 touching (●),
then states "safe to run fully in parallel: WS-1 + WS-3 + WS-7 + WS-8". Both
cannot hold — WS-3 needs fold-safe evaluation plumbing from that file. Either
WS-3 confines itself to a new `bench/hybrid.py` with read-only imports from
`harness.py`, or the two are sequenced.

### 7.6 🟢 Strengthens §3.1 — the embedding is *structurally* blind to capsule

§3.1 treats the reviewer's point as a hypothesis to test. It is stronger than
that. `pf-bacteria.genome_vector` is MDS over **PanACoTA core-genome distances**
(`current_status.md` §6.1). Capsule and O-antigen loci are hypervariable and
horizontally transferred — routinely *excluded* from core-genome alignments by
construction. The embedding is not merely insensitive to capsule; it is built
from a signal that omits it. This predicts WS-1 should work, and is worth saying
in any write-up.

### 7.7 🟠 Amends WS-1 — statistical guards on the two proposed features

**`defence_mismatch` needs a prevalence band.** `MazEF` is in 338/403 strains
(84%) and `RM_Type_I` in 327 (81%). Depletion is unmeasurable for systems that
near-universal, so a raw fraction will be noise. Whatever signal exists lives in
the rarer systems — apply a prevalence filter and a significance test, not a
fraction. (The `significant_terms` machinery Stage C already uses does exactly
this job.)

**`receptor_compat` is `breadth` stratified by receptor match.** That is a good
feature and a fair test of the reviewer's hypothesis, but it is not biophysical
compatibility, and it will correlate with `breadth` (weighted 0.25) and with the
`generalist` baseline. Check collinearity when fitting, or the conclusion
"receptor biology matters" may be partly a re-weighting of host-range breadth.

### 7.8 ⚪ Amends §3.5 — cite the prior-art figures or drop the digits

PhageHostLearn "AUC 81.8%", HostPhinder "81% genus", "a 2024 benchmark reviewed
**27** such tools". The tools are real and the positioning is right, but these
specific figures are unverified here and are exactly what a judge spot-checks.
Attach citations or state them qualitatively.

### 7.9 Revised order for §6

1. **Stage D plumbing** (§7.4) — ~1 h. Without it WS-1 measures nothing.
2. **Per-strain scores + paired significance** (§7.3) — ~1 h. Without it WS-1's
   result cannot be read.
3. **WS-1.**
4. Then WS-3 / WS-4 / WS-7 as capacity allows.

Steps 1 and 2 are prerequisites, not polish: they are the difference between WS-1
producing a finding and WS-1 producing another 0.633.

### 7.10 ✅ §7.3 landed — and it settles §2.1, §2.2 and §3.1/§3.2

`bench/significance.py` now reports a paired bootstrap CI and a Wilcoxon p-value
per method, paired by strain, and `benchmark.json` carries the per-strain vectors
(`per_strain`) so any comparison can be recomputed without a rerun. Measured over
390 strains, `infect_threshold = 0.0`:

**vs `funnel`** — does the stage carry signal?

| Method | ΔP@10 | 95% CI | p | W/L/T |
|---|--:|:--|--:|:--|
| `+ receptor & defence` | **+0.0136** | +0.0056 … +0.0218 | 0.002 | 92/56/242 |
| `+ receptor only` | **+0.0131** | +0.0062 … +0.0200 | 7.9e-05 | 77/43/270 |
| `- RBP arm` | +0.0010 | −0.0031 … +0.0051 | 0.71 | 30/27/333 |
| `+ defence only` | +0.0008 | −0.0056 … +0.0072 | 0.82 | 61/57/272 |
| `- neighbour prior` | −0.0997 | −0.1172 … −0.0831 | 5.6e-25 | 57/212/121 |

**vs `phylo_nn`** — do we beat the baseline that matters? `funnel` **+0.0479**
(CI +0.0269 … +0.0692, p = 1.5e-06).

Four conclusions, each now backed by an interval rather than a decimal place:

1. **§3.1 is confirmed, and it is the receptor half.** Promoting receptor
   compatibility into scoring buys +0.0136 P@10 with the CI clear of zero. The
   reviewer identified a real defect. It is a modest effect, not a
   transformation — quote the interval, not the point estimate.
2. **§3.2's defence half does nothing detectable** (+0.0008, CI spans zero,
   p = 0.82). Consistent with §7.7's warning: the common systems are too
   near-universal to carry signal. Worth reporting as a negative result, and
   worth retrying only with a prevalence filter.
3. **§2.2's conclusion survives, now honestly.** With the RBP arm genuinely
   wired (§7.1/§7.4), `- RBP arm` is +0.0010 with CI −0.0031 … +0.0051. The arm
   contributes nothing **detectable** — a measured null, not an artefact of
   disconnected plumbing. §2.2's reasoning was wrong; its answer was right.
4. **§2.1's alarm can be stood down, with a caveat.** The +4.8 points over
   `phylo_nn` is real (p = 1.5e-06), but the CI runs from +2.7 to +6.9 points.
   Quote the range. Separately, `generalist` vs `phylo_nn` is −0.0041, CI
   spanning zero — **the two baselines are statistically indistinguishable**, so
   "beats three baselines" is really "beats two distinct ones".

Cost: the benchmark is now 174 s rather than 69 s, mostly from the three added
funnel configurations rather than from the bootstrap.

**Consequence for §6:** WS-1 has reported. Its receptor half should ship; its
defence half should not, pending a prevalence filter. Anything built on Stage D
from here should be measured with a CI, not a mean.
