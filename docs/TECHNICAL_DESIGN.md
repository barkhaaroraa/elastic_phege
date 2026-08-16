# PhageForge — Technical Design

> **Status:** design document, pre-implementation
> **Audience:** build team (Python, Elasticsearch, ML)
> **Companion docs:** [`ARCHITECTURE.md`](./ARCHITECTURE.md) (system design, mappings, query bodies) · [`SOLUTION_OVERVIEW.md`](./SOLUTION_OVERVIEW.md) (non-technical review)

This document covers **what we are building, why, how we will know it works, and what will go wrong**. The system design itself lives in `ARCHITECTURE.md` and is referenced rather than duplicated.

---

## 1. Context

Built for Elastic's **Forge the Future 2026** hackathon (Elastic × AWS, hosted on HackCulture; theme *"Ground it. Automate it. Ship it."*; INR 600,000 prize pool; live stage at Elastic{ON} Mumbai, 30 Sept). Elastic must be the primary backend.

Relevant tracks:
- **Competitive Benchmarking** — primary. Our benchmark table is the main deliverable.
- **Intelligent Search and AI Platforms** — the retrieval funnel itself.
- **Observability and Autonomous Operations** — secondary, via APM instrumentation.
- **Industry Solutions and Vertical Experiences** — life sciences framing.

> ⚠️ **Open item.** Public listings state a submission deadline of **13 Aug 2026**. Confirm the real date against your registration confirmation before committing to §9's schedule; the milestone table assumes a runway and must be compressed or rebased if the deadline has passed.

---

## 2. Goals and non-goals

### Goals

| # | Goal | Measured by |
|---|---|---|
| G1 | Reduce a ~10,000-phage search space to a ranked top-10 for an uncharacterised strain | Funnel returns ≤ 10 with reasons, p95 < 2 s for a known strain |
| G2 | Beat naive baselines on held-out strains | precision@10, recall@10 vs. 3 baselines (§7) |
| G3 | Every prediction carries a **biological reason** and an **evidence tier** | 100% of returned phages have both populated |
| G4 | Make Elastic load-bearing at every funnel stage | No scoring logic in Python except set-cover |
| G5 | Design a receptor-diverse cocktail, not just a top-10 list | Cocktail endpoint + rationale |
| G6 | Accept a novel genome (FASTA) end to end | Upload → features → prediction without manual steps |

### Non-goals — state these out loud in the demo

- **Not** beating Gaborieau or PhageHostLearn on raw AUROC. They are Nature papers with dedicated models; we are the decision-support and prioritisation layer that does not exist around them.
- **Not** a clinical treatment recommender. No dosing, no *in vivo* efficacy, no patient-level advice.
- **Not** a general-purpose phage annotation pipeline (Sphae, PhageScope already do that well).
- **Not** claiming strain-level rigour outside *Escherichia* and *Klebsiella*.

### The positioning that wins

Do not sell *"our AI predicts phage–host interactions."* You will be benchmarked against two Nature papers and lose in three weeks.

Sell **the decision-support layer**: one queryable substrate unifying genomes, interaction matrices, biological features and literature — with provenance, explainability, and experiment planning on top. Nobody has built that, and it is exactly what a phage biobank needs.

---

## 3. Stack

| Layer | Choice | Rationale |
|---|---|---|
| Backend | **Python 3.11 + FastAPI** | Bioinformatics tooling is Python; thin orchestration only |
| Data layer | **Elastic Cloud** (9.x) | ELSER, Agent Builder, LTR, BBQ available without self-hosting |
| Vectors | ESM-2 650M (1280-dim), MinHash k-mer sketches (256-dim) | §5 |
| Ranker | XGBRanker → Eland LTR export → in-cluster rescorer | Native Elastic inference |
| NL interface | Elastic Agent Builder (ES\|QL tools) | GA; second surface on the same data |
| Frontend | Next.js (App Router) + Tailwind | Fast to build, good for the funnel viz |
| Observability | Elastic APM → Kibana | Secondary track, ~1 day of work |
| Orchestration | Prefect or plain Makefile + scripts | Ingest is batch, not streaming |

---

## 4. Data sources and ingest

### 4.1 Per-source plan

| Source | Records (expected) | Format | Licence | Tier | Extraction |
|---|---|---|---|:--:|---|
| **Zenodo 13831957** (Gaborieau) | 403 strains × 96 phages = 38,688 interactions, **with negatives** | tar.gz, 146.8 MB | CC (verify) | 1 | Parse matrix CSV → `pf-interactions`; genomes → `pf-bacteria`/`pf-phages` |
| **PhageHostLearn** + Zenodo 11061100 | Klebsiella matrix, RBPs, K-loci | notebooks + processed data | MIT (code) | 1 | Reuse their extraction notebooks; do not reimplement |
| **INPHARED** | ~30k phage genomes + taxonomy TSV + MASH sketches | FASTA, TSV | AGPL-3.0 | — | Bulk → `pf-phages` (exploratory tier) |
| **PhageScope / MVP / Virus-Host DB** | bulk interactions | TSV/API | varies | 2 | Normalise to `pf-interactions` with `resolution: species` |
| **NCBI / BV-BRC** | bacterial assemblies | FASTA | public | — | Genomes for `pf-bacteria` |
| **Europe PMC OA** | full text | REST/XML | OA subset | 3 | Chunk → `pf-literature` via `semantic_text` |

### 4.2 The three data caveats — design around them, do not paper over them

1. **Public interactions outside the two matrices are positives-only and species-level.** The imagined "100,000+ strain-level interactions" is really **~40k good ones** plus a large pile of lower-resolution evidence. Hence `resolution` as a first-class field and the evidence-tier system.
2. **No negatives means no ranking.** A ranker must know what *doesn't* infect. Only Tier 1 has that, which is why training and evaluation are confined to Tier 1 and why organism scope is two genera.
3. **ESM-2 embedding is compute-heavy.** Precompute offline, index vectors, **never embed at query time** except for a single uploaded genome's receptor proteins.

### 4.3 Ingest mechanics

- `_bulk` in 5,000-doc batches, `refresh=false` during load, single explicit refresh at the end.
- Index templates created before load; no dynamic mapping (a dynamically-mapped `dense_vector` will not be indexed correctly).
- Idempotent: `_id` derived deterministically (`sha1(source + native_id)`) so re-runs upsert.
- Every ingested doc records `source` and `ingested_at` for provenance and reproducible benchmark snapshots.

---

## 5. Feature extraction pipeline

| Target | Tool | Output field(s) | Notes |
|---|---|---|---|
| Bacterial K-locus / capsule | **Kaptive 2.x** | `k_locus`, `k_locus_conf` | Klebsiella-first; keep the confidence call |
| Bacterial MLST | mlst / PubMLST scheme | `st` | ST11, ST258 etc. |
| O-antigen / LPS loci | Kaptive O-typing + gene screen | `o_antigen`, `lps_genes` | Primary receptor class for *E. coli* |
| Antiphage defence | **DefenseFinder** or PADLOC | `defence_systems` | Marginal predictive value — kept for explanation, not lift |
| AMR genes | **AMRFinderPlus** | `amr_genes` | Context for the researcher; not a model feature |
| Phage lifestyle | **BACPHLIP** | `lifestyle`, `lifestyle_conf` | Temperate → down-rank, not drop |
| Phage RBPs | **PhageRBPdetect 2.x** | `pf-proteins` role=RBP | The single most important feature |
| Depolymerases | HMM screen | `has_depolymerase` | Capsule-degrading, matters for Klebsiella |
| Protein embeddings | **ESM-2 650M** | `esm2_vector` (1280) | See risk R1 |
| Genome sketches | MinHash / k-mer profile | `genome_vector` (256) | Cheap tier + fallback |

**Receptor-class assignment** is derived, not measured: map RBP homology and annotation keywords onto the controlled vocabulary (`LPS | K-capsule | OmpC | OmpF | BtuB | LamB | TonB | flagellar | pili | unknown`). This is approximate and must be labelled as such — it drives cocktail diversity, so a reviewer should be asked whether the vocabulary is right (see `SOLUTION_OVERVIEW.md` §12).

---

## 6. Elastic data layer

Full mappings and query bodies are in [`ARCHITECTURE.md` §5 and §7](./ARCHITECTURE.md). Implementation notes that belong here:

- **Verified available:** `rrf` retriever · `learning_to_rank` rescorer (native since 8.12; Eland ships an XGBRanker export helper) · `significant_terms` · `semantic_text` + ELSER · `dense_vector` up to 4096 dims · `int8_hnsw` and `bbq_hnsw` quantization.
- **`learning_to_rank` is technical preview.** Pin the Elastic Cloud version. Fallback if it misbehaves: run XGBoost in FastAPI over the 500 survivors and label the architecture honestly — but this weakens the "Elastic is load-bearing" claim, so treat it as a last resort (risk R5).
- **`window_size` must be ≥ `from + size`** on the rescorer, and should equal the candidate-set size (500) so every survivor is scored.
- **Quantization assignment:** `bbq_hnsw` on the 1280-dim ESM-2 vectors (BBQ performs best at high dimensionality, ~32× compression); `int8_hnsw` on the 256-dim sketches. Use `rescore_vector` with full-precision oversampling where recall matters.
- **Denormalise the phage's primary RBP embedding onto `pf-phages`** as `rbp_match_vector` so Stage A is one round trip instead of a protein→phage join.

---

## 7. Model training and evaluation

**This section is the primary deliverable.** Given judging on "real performance, reliability, and impact, not pitch polish," a defensible benchmark beats any amount of UI.

### 7.1 Training data

Tier 1 only: Gaborieau (403 × 96, with negatives) + PhageHostLearn Klebsiella. Positives = `infects: true`; negatives = explicitly measured non-infection. Species-level records are **excluded from training** — they would inject label noise at exactly the resolution we claim to operate at.

### 7.2 Splits — the thing most likely to be got wrong

**Split by strain, never by interaction.** A strain must appear in exactly one fold. Splitting on the 38,688 pairs randomly leaks a strain's other 95 measurements into training and produces a beautiful, meaningless AUROC.

- 5-fold grouped CV, groups = `strain_id`.
- A held-out set of strains never touched during development, used once for the final number.
- Report per-organism results separately (*Escherichia* vs. *Klebsiella*) — cross-organism generalisation is a claim, so it needs its own number.

### 7.3 Model

`XGBRanker` (`rank:pairwise`), groups = strain, features per §7.4 of `ARCHITECTURE.md`. Exported via Eland's LTR helper to `phageforge-ltr-v1`. Keep the model small — inference runs inside the cluster per rescore window.

### 7.4 Baselines — all three, no exceptions

| Baseline | Why it matters |
|---|---|
| **Random** | Floor. |
| **Most-generalist phage** | The one that actually competes. Ranking every strain by global phage host-range breadth is a strong, dumb baseline; if we do not clearly beat it, the model is learning nothing strain-specific. |
| **Phylogenetic nearest neighbour** | Take the single most similar tested strain, copy its results. The "obvious" method a microbiologist would use by hand. |

### 7.5 Metrics

- **precision@10** — of our 10 recommendations, how many are true hits? *The headline number.*
- **recall@10** — of all phages that truly infect this strain, what fraction did we surface in 10?
- **AUROC / AUPRC** — comparability with published work (AUPRC matters more; positives are rare).
- **Enrichment over random** — the demo slide: *"random gets 1 in 10, we get 7 in 10."*
- **Cocktail coverage** — fraction of held-out strains covered by a 4-phage cocktail vs. the top-4 individually ranked phages.

### 7.6 The benchmark table (target output)

| Method | P@10 | R@10 | AUPRC | Cocktail coverage |
|---|:--:|:--:|:--:|:--:|
| Random | | | | |
| Most-generalist phage | | | | |
| Phylogenetic NN | | | | |
| **PhageForge funnel** | | | | |
| PhageForge − RBP features (ablation) | | | | |
| PhageForge − neighbour prior (ablation) | | | | |

Ablations are included deliberately: they show *which part* of the funnel carries the signal, which is a far more credible claim than a single aggregate number.

---

## 8. Implementation notes by component

### 8.1 Funnel orchestration
One module, `phageforge/funnel.py`, with a function per stage returning `(results, stage_metadata)`. Stage metadata (counts, latency, query used) is what gets persisted to `pf-predictions` and rendered in the UI. Never let a stage silently return zero — an empty stage must surface as an explicit "no candidates passed filter X" state, not a blank shortlist.

### 8.2 Genome upload path
FASTA → QC (contig count, N50, completeness) → reject with a clear message if below threshold → Kaptive/MLST/DefenseFinder → ESM-2 on receptor proteins → ephemeral `pf-bacteria` doc → funnel. Runs as a background job with `/jobs/{id}` polling; expect minutes, not seconds. **Do not block the HTTP request on this.**

### 8.3 Cocktail
Greedy set-cover per `ARCHITECTURE.md` §7.6. Keep λ configurable and expose it in the UI — a reviewer will want to see how the set changes as receptor-diversity pressure increases.

### 8.4 Plate layout export
96-well CSV: shortlist phages × dilution series, positive/negative controls included, in the layout a bench scientist would actually pipette. Small feature, disproportionate credibility — it shows the tool understands the downstream experiment. Ask the reviewer to sanity-check the layout.

### 8.5 Frontend
The **funnel visualisation is the demo**. Show 10,000 → 500 → 50 → 10 as counts narrowing with each stage clickable to reveal what it did and why. Evidence-tier badges everywhere. Resist building a generic dashboard.

---

## 9. Milestones

Six working weeks, assuming a runway from mid-August. **Rebase against the confirmed deadline (§1) before use.** If the runway is shorter, cut in this order: Agent Builder → cocktail UI → second organism → observability. Never cut the benchmark.

| Week | Deliverable | Done when |
|:--:|---|---|
| **W1** | Ingest + indices | All six indices populated from Gaborieau + PhageHostLearn + INPHARED; counts recorded |
| **W2** | Feature extraction + embeddings | Kaptive/BACPHLIP/RBP outputs indexed; ESM-2 vectors in `pf-proteins`; sketch fallback proven |
| **W3** | Funnel stages A–C | RRF candidate generation, neighbour prior, `significant_terms` explanations returning sane biology |
| **W4** | LTR + **benchmark** | XGBRanker trained on strain-grouped folds, deployed via Eland, benchmark table populated against all 3 baselines |
| **W5** | Cocktail + UI + Agent Builder | Set-cover, funnel viz, plate export, NL query surface |
| **W6** | Hardening, APM, demo, writeup | APM spans + Kibana dashboard; end-to-end demo rehearsed; README + architecture published |

**W4 is the critical path.** If the benchmark is not in hand by end of W4, cut scope elsewhere immediately.

---

## 10. Risks and mitigations

| # | Risk | Impact | Mitigation |
|:--:|---|---|---|
| **R1** | **ESM-2 embedding compute** on day one | High — blocks W2 | Budget GPU early; MinHash/k-mer sketches are a proven fallback that index into Elastic identically and are ~100× cheaper. Decide by end of W2, do not drift. |
| **R2** | **Sparse negatives** outside Tier 1 | High — no ranking signal | Confine training/eval to Tier 1; be explicit in the docs and demo. Do not synthesise negatives by assuming untested = non-infecting. |
| **R3** | **Elastic becomes a metadata store** | Fatal for judging | Architectural rule: no scoring in Python except set-cover. Review every stage against this. |
| **R4** | **Elastic Cloud trial limits** (ELSER + ESM-2 + LTR concurrently) | Medium | Size the deployment early; ELSER inference at ingest only; consider trimming the literature corpus before the ML node budget bites. |
| **R5** | **`learning_to_rank` is technical preview** | Medium | Pin the version; fallback is Python XGBoost over 500 survivors (weakens R3 — last resort). |
| **R6** | **Scope creep to more organisms** | Medium | Two organisms rigorously beats five vaguely. The scope boundary in `ARCHITECTURE.md` §9 is a commitment. |
| **R7** | **Receptor-class vocabulary is approximate** | Medium — drives cocktail logic | Get it reviewed by a microbiologist (that is what `SOLUTION_OVERVIEW.md` §12 is for); label as heuristic in the UI. |
| **R8** | **Deadline already passed** (§1) | Existential for the hackathon | Confirm immediately; the system still stands as a portfolio/product artifact regardless. |

---

## 11. Definition of done

- [ ] Six indices populated, record counts documented as actuals (not estimates)
- [ ] Funnel returns ≤ 10 ranked phages for any *Escherichia* or *Klebsiella* strain, p95 < 2 s
- [ ] Every result carries a biological reason string and an evidence tier
- [ ] Benchmark table populated with strain-grouped held-out splits and all 3 baselines + 2 ablations
- [ ] Cocktail endpoint returns a receptor-diverse set with stated rationale
- [ ] FASTA upload → prediction works end to end without manual intervention
- [ ] Agent Builder answers the four example NL queries
- [ ] APM spans per stage visible in a Kibana dashboard
- [ ] `SOLUTION_OVERVIEW.md` reviewed by a domain expert, feedback incorporated
- [ ] Demo rehearsed end to end in under 5 minutes
