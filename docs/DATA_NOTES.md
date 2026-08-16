# Data Notes — what is actually in the Gaborieau archive

> Recorded at ingest time from the extracted archive, not from the paper.
> Where reality differs from `ARCHITECTURE.md`, **this file is correct** and the
> architecture doc is aspirational.

**Source:** Zenodo record [13831957](https://zenodo.org/records/13831957),
`coli_phage_interactions_2023.tar.gz`
**Verified:** MD5 `380d288b8cd7629660603447b9856fc9` ✅ · 146,786,594 bytes ✅ · gzip OK ✅
**Safety audit:** 862 entries (803 files, 59 dirs); zero absolute paths, zero `../`
traversal, zero symlinks/devices, zero setuid. Extracted with `tarfile` `filter="data"`.

---

## What we got

| Asset | Location | Status |
|---|---|:--:|
| Interaction matrix | `data/interactions/interaction_matrix.csv` | ✅ |
| Bacterial strain features | `data/genomics/bacteria/picard_collection.csv` | ✅ |
| Defence systems (per strain) | `data/genomics/bacteria/defense_finder/` | ✅ |
| O-type / serotype | `data/genomics/bacteria/isolation_strains/o_type/output.tsv` | ✅ |
| LPS type | `data/genomics/bacteria/.../outer_core_lps/LPS_type_waaL_370.txt` | ✅ |
| Outer-membrane proteins | `data/genomics/bacteria/.../outer_membrane_proteins/` | ✅ |
| Capsules (Kaptive output) | `data/genomics/bacteria/capsules/` | ✅ |
| Strain phylogenetic distances | `.../panacota/tree/370+host_distance_matrix_orignames.tsv` | ✅ |
| Phage metadata | `data/genomics/phages/guelin_collection.csv` | ✅ |
| Phage genomes | `data/genomics/phages/FNA/*.fna` (97 files) | ✅ |
| Phage RBP identifiers | `data/genomics/phages/RBP/RBP_list.csv` | ⚠️ IDs only |
| **Bacterial genomes** | — | ❌ **absent** |
| **Protein sequences (any)** | — | ❌ **absent** |

The 33 `.nsq/.nin/.nhr` files under `isolation_strains/Phylogroup/db/` are BLAST
databases, not retrievable sequences.

---

## The interaction matrix

Semicolon-delimited. Header is `bacteria;<phage_1>;…;<phage_96>`.

- **402 strains × 96 phages = 38,592 cells**
- 157 cells empty (untested) → **38,435 measured**
- The paper reports 403 × 96 = 38,688. `picard_collection.csv` lists 403 strains,
  so **one strain has metadata but no interaction row**. Minor, but our counts are
  38,435 measured, not 38,688 — report ours, not the paper's.

### Values are graded 0–4, not binary

| Score | Cells | Share |
|:--:|--:|--:|
| 0.0 | 30,459 | 79.2% |
| 1.0 | 3,014 | 7.8% |
| 2.0 | 2,306 | 6.0% |
| 3.0 | 1,368 | 3.6% |
| 4.0 | 1,288 | 3.4% |
| *(empty)* | 157 | — |

This is a **graded lysis score**, not EOP. We store the raw value in
`score_raw` and derive `infects` with a configurable threshold.

**Default threshold: `score_raw > 0`.** Grounded in the authors' own analysis
code, where `>0` / `>0.0` comparisons on the matrix dominate by a wide margin.
That yields **7,976 positives (20.7%)** and 30,459 negatives.

> ⚠️ This is precisely review question 3 in `SOLUTION_OVERVIEW.md` §12. A `1.0`
> may be turbid/faint lysis — possibly lysis-from-without rather than productive
> infection. If our reviewer says the cut belongs at ≥2, that single change moves
> positives from 7,976 to 4,962 and will move every benchmark number. The
> threshold is a config value (`PF_INFECT_THRESHOLD`) for exactly this reason.

**Negatives are real.** 30,459 measured non-infections is what makes ranking
possible and is why this dataset is the whole Tier-1 foundation.

---

## Consequences for the architecture

Three deviations from `ARCHITECTURE.md`, all forced by what is actually shipped.

### 1. Bacterial `genome_vector` — from phylogeny, not k-mers

No bacterial genomes, so no k-mer sketch. Instead we derive strain embeddings
from `370+host_distance_matrix_orignames.tsv` (404 strains, symmetric, PanACoTA
distances) by classical MDS into `SKETCH_DIMS`.

This is arguably *better* than the original plan: kNN over these vectors
reproduces genuine phylogenetic neighbourhood, which is exactly what Stage B's
neighbour-transfer needs — rather than approximating it with k-mer similarity.

### 2. Phage `genome_vector` — as designed

MinHash/k-mer sketches over the 97 phage FASTAs. Single-contig complete genomes
(`>AN17_P8`, one record each), so this is clean.

### 3. RBP embeddings — must gene-call first

`RBP_list.csv` gives 137 rows / **130 RBP identifiers across 89 phages**, typed as
`fiber` (104), `spike` (23), `confirmed` (3), `NA` (7). Identifiers are locus tags
in Prokka form — `AN17_P8_00026` = 26th CDS of phage `AN17_P8`.

Sequences are **not** shipped, so the RBP arm requires gene-calling the phage
FASTAs (`pyrodigal` 3.7.1, meta mode) and recovering the CDS by ordinal position.

### Verification result (measured, 89 phages / 130 RBP ids)

The assumption is **broadly correct but not exact.**

Strong supporting evidence: RBP-indexed proteins have **median 683 aa** against
**median 128 aa** for all 8,536 called CDS — a 5.3× enrichment, precisely the
size profile of tail fibres and spikes. The ordinals are clearly pointing at the
right *kind* of protein.

But the mapping is not one-to-one:

| Outcome | Count | Share |
|---|--:|--:|
| Ordinal lands on a plausible RBP (≥ 300 aa) | 111 | 85.4% |
| Ordinal lands on something too short | 16 | 12.3% |
| Ordinal exceeds the CDS count entirely | 3 | 2.3% |
| *(of the 16 short ones, recoverable within ±3)* | *15* | — |

The offsets run in **both directions** (`409_P1_00051` → CDS 49; `409_P3_00008` →
CDS 9), so this is not a single systematic shift. The likely cause is that Prokka
numbers tRNAs and other non-CDS features in the same `locus_tag` sequence, and
its Prodigal settings differ from pyrodigal's meta mode — so the two gene calls
drift apart at different points in each genome.

### Decision: confidence-tiered, never silently guessed

A ±3 window recovers 126/130 (96.9%), but "longest protein nearby" could just as
easily grab a portal or terminase as a fibre — a silent wrong-protein error,
which is worse than a missing one. So RBP proteins carry an explicit confidence:

* **`high`** — exact ordinal, protein ≥ 300 aa. **111 RBPs.** Safe to embed.
* **`heuristic`** — exact ordinal implausible, recovered within ±3. **15 RBPs.**
  Embedded but flagged; Stage A may down-weight or exclude them.
* **`unresolved`** — ordinal out of range. **3 RBPs.** Not embedded.

This mirrors the evidence-tier principle used everywhere else in the project: the
uncertainty is recorded on the document, not averaged away.

> The clean fix is to run Prokka itself and reproduce the authors' exact
> numbering. That is a heavy external dependency and is deferred; the tiering
> above is the honest interim.

### 4. No bacterial receptor vectors — Stage A's RBP arm changes shape

There are no bacterial protein sequences, so the designed
"phage RBP vector ↔ bacterial receptor vector" match is **not computable**.

Bacterial surface biology is available only as *categorical* features (O-type,
LPS type, capsule type, OMP presence/absence). Those remain excellent for Stage C
`significant_terms` and as ranker features — but they are not vectors.

**Adapted Stage A:** the third retriever arm becomes RBP-embedding similarity
between candidate phages and the RBPs of phages already known to infect strains
near the target. That is still a genuine vector arm and is well-founded — it is
collaborative filtering in RBP space rather than a direct biophysical match.

---

## Strain metadata columns (`picard_collection.csv`, `;`-delimited)

`bacteria` · `Gembase` · `Host` · `Origin` · `Pathotype` · `Clermont_Phylo` ·
`ST_Warwick` · `O-type` · `H-type` · `Mouse_killed_10` · `Capsule_ABC` ·
`Capsule_GroupIV_e` · `Capsule_GroupIV_e_stricte` · `Capsule_GroupIV_s` ·
`Capsule_Wzy_stricte` · `LPS_type` · `Collection` · `Klebs_capsule_type` ·
`n_defense_systems` · `n_infections` · `ABC_serotype`

This single file covers most of `pf-bacteria`; defence-system *names* come from
`defense_finder/370+host_defense_systems_subtypes.csv` (wide, `;`-delimited,
one column per system, counts — convert to a list where count > 0).

## Phage metadata columns (`guelin_collection.csv`, `;`-delimited)

`phage` · `Morphotype` · `Family` · `Genus` · `Species` · `Genome_size` ·
`Phage_host` · `Phage_host_phylo` · `Old_Family` · `Subfamily` · `Old_Genus`

96 phages. Note there is **no lifestyle (virulent/temperate) column** — BACPHLIP
is not run here, so `lifestyle` stays `unknown` until we add it.
