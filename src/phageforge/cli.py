"""Command-line entry points. ``python -m phageforge.cli <command>``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from phageforge import config
from phageforge.es import bulk, mappings


def _client():
    client = config.get_client()
    health = config.require_healthy(client)
    print(f"cluster: {health['cluster_name']} [{health['status']}]")
    return client


def cmd_indices(args: argparse.Namespace) -> int:
    client = _client()
    actions = mappings.create_indices(client, recreate=args.recreate)
    for name, action in actions.items():
        print(f"  {name:20s} {action}")
    return 0


def cmd_ingest(args: argparse.Namespace) -> int:
    from phageforge.ingest import gaborieau, zenodo

    fetched = zenodo.download()
    if not fetched.verified:
        print(
            f"ERROR: archive failed verification (md5={fetched.md5}, size={fetched.size})",
            file=sys.stderr,
        )
        return 1
    print(f"archive verified: {fetched.size:,} bytes, md5 {fetched.md5}")
    zenodo.extract()

    client = _client()
    mappings.create_indices(client, recreate=args.recreate)

    interactions, stats = gaborieau.read_interactions()
    print(
        f"\ninteraction matrix: {stats.strains} strains x {stats.phages} phages "
        f"= {stats.cells:,} cells"
    )
    print(f"  measured : {stats.measured:,}   untested (dropped): {stats.untested}")
    print(
        f"  positives: {stats.positives:,} ({stats.positives / stats.measured:.1%})"
        f"   negatives: {stats.negatives:,}"
    )
    print(f"  threshold: score > {gaborieau.INFECT_THRESHOLD}")

    bacteria = gaborieau.read_bacteria()
    phages = gaborieau.read_phages(interactions)

    plan = [
        (
            config.INTERACTIONS,
            (
                (bulk.doc_id(gaborieau.SOURCE, f"{d['phage_id']}|{d['host_id']}"), d)
                for d in interactions
            ),
        ),
        (
            config.BACTERIA,
            ((bulk.doc_id(gaborieau.SOURCE, d["strain_id"]), d) for d in bacteria),
        ),
        (
            config.PHAGES,
            ((bulk.doc_id(gaborieau.SOURCE, d["phage_id"]), d) for d in phages),
        ),
    ]

    print()
    for index, docs in plan:
        n = bulk.bulk_index(client, index, docs)
        print(f"  {index:20s} indexed {n:,}  ->  {bulk.count(client, index):,} in index")

    print("\ncounts by evidence tier / resolution:")
    for field in ("evidence_tier", "resolution"):
        agg = client.search(
            index=config.INTERACTIONS,
            size=0,
            aggs={"by": {"terms": {"field": field, "size": 10}}},
        )
        buckets = agg["aggregations"]["by"]["buckets"]
        print(f"  {field:15s} " + ", ".join(f"{b['key']}={b['doc_count']:,}" for b in buckets))
    return 0


def cmd_features(args: argparse.Namespace) -> int:
    """Compute and attach genome vectors (Phase 3 items 1 and 2)."""
    import numpy as np

    from phageforge.features import phylo, sketch
    from phageforge.ingest import gaborieau

    client = _client()

    # --- phages: tetranucleotide-frequency sketches -------------------------
    print("\nphage TNF sketches")
    vectors = sketch.sketch_directory(gaborieau.FNA_DIR)
    indexed = {
        h["_source"]["phage_id"]: h["_id"]
        for h in _scan_ids(client, config.PHAGES, "phage_id")
    }
    matched = {pid: vec for pid, vec in vectors.items() if pid in indexed}
    print(f"  sketched {len(vectors)} FASTA files; {len(matched)} match indexed phages")
    missing = sorted(set(indexed) - set(matched))
    if missing:
        print(f"  WARNING: {len(missing)} indexed phages have no FASTA: {missing[:5]}")
    applied = bulk.bulk_update(
        client,
        config.PHAGES,
        ((indexed[pid], {"genome_vector": vec.tolist()}) for pid, vec in matched.items()),
    )
    print(f"  attached genome_vector to {applied} phages")

    # --- bacteria: classical MDS over the PanACoTA distance matrix ----------
    print("\nbacteria phylogenetic embeddings (classical MDS)")
    ids, matrix = gaborieau.read_distance_matrix()
    result = phylo.classical_mds(np.array(matrix), ids, config.SKETCH_DIMS)
    quality = phylo.reconstruction_error(np.array(matrix), result.coords, result.scale)
    print(f"  matrix {len(ids)}x{len(ids)} -> {result.coords.shape}")
    print(f"  coordinate scale   : x{result.scale:.1f} (order-preserving)")
    print(f"  variance explained : {result.explained:.1%}")
    print(f"  negative eigenmass : {result.negative_mass:.2%}")
    print(
        f"  reconstruction     : r={quality['pearson_r']:.4f} "
        f"relative MAE={quality['relative_mae']:.2%}"
    )
    if quality["pearson_r"] < 0.95:
        print("  ERROR: embedding does not reproduce the phylogeny", file=sys.stderr)
        return 1

    strain_ids = {
        h["_source"]["strain_id"]: h["_id"]
        for h in _scan_ids(client, config.BACTERIA, "strain_id")
    }
    pairs = [
        (strain_ids[sid], {"genome_vector": result.coords[i].tolist()})
        for i, sid in enumerate(result.ids)
        if sid in strain_ids
    ]
    print(f"  {len(pairs)} of {len(result.ids)} matrix strains match indexed strains")
    unmatched = sorted(set(strain_ids) - set(result.ids))
    if unmatched:
        print(f"  WARNING: {len(unmatched)} indexed strains absent from matrix: {unmatched[:5]}")
    applied = bulk.bulk_update(client, config.BACTERIA, iter(pairs))
    print(f"  attached genome_vector to {applied} strains")

    print("\ncoverage:")
    for index in (config.PHAGES, config.BACTERIA):
        with_vec = bulk.count(client, index, {"exists": {"field": "genome_vector"}})
        print(f"  {index:20s} {with_vec:,} / {bulk.count(client, index):,} have genome_vector")
    return 0


def cmd_proteins(args: argparse.Namespace) -> int:
    """Recover RBPs by gene calling and embed them with ESM-2 (Phase 3 item 3)."""
    from phageforge.features import proteins
    from phageforge.ingest import gaborieau

    client = _client()
    # Only pf-proteins is ever recreated here: it is derived entirely from this
    # command, whereas pf-phages holds ingested data this command only annotates.
    mappings.create_indices(client, (config.PROTEINS,), recreate=args.recreate)
    mappings.create_indices(client, (config.PHAGES,))

    rbp_index = gaborieau._rbp_index()
    n_ids = sum(len(v) for v in rbp_index.values())
    print(f"\nRBP_list.csv: {n_ids} identifiers across {len(rbp_index)} phages")

    print("gene calling (pyrodigal, meta mode)")
    records, stats = proteins.resolve_rbps(rbp_index, gaborieau.FNA_DIR)
    print(f"  called {stats.cds_called:,} CDS across {stats.phages_called} phages")
    print(
        f"  median length: RBP-indexed {stats.median_rbp_aa():.0f} aa "
        f"vs {stats.median_cds_aa():.0f} aa for all CDS "
        f"({stats.median_rbp_aa() / max(stats.median_cds_aa(), 1):.1f}x)"
    )
    print("  confidence tiers:")
    for tier in (proteins.HIGH, proteins.HEURISTIC, proteins.UNRESOLVED):
        n = stats.by_confidence.get(tier, 0)
        print(f"    {tier:12s} {n:4d}  ({n / max(stats.total, 1):.1%})")
    for record in records:
        if record.confidence == proteins.UNRESOLVED:
            print(f"      skipped {record.rbp_id}: {record.note}")

    embeddable = [r for r in records if r.gene is not None]
    if not embeddable:
        print("ERROR: nothing resolved to embed", file=sys.stderr)
        return 1
    if args.resolve_only:
        print("\n--resolve-only: stopping before embedding")
        return 0

    print(f"\nESM-2 embeddings ({config.ESM_MODEL}, {config.ESM_DIMS} dims, CPU)")
    sequences = {r.protein_id: proteins.clean_sequence(r.gene.protein) for r in embeddable}
    raw = proteins.embed_proteins(sequences, batch_size=args.batch_size)
    print(f"  embedded {len(raw)} proteins")

    vectors, mean = proteins.center_vectors(raw)
    path = proteins.save_mean_vector(mean)
    spread = _cosine_spread(vectors)
    print(
        f"  centred against the corpus mean (saved to {path.name}); "
        f"pairwise cosine now {spread[0]:.2f}..{spread[1]:.2f}, median {spread[2]:.2f}"
    )

    now = gaborieau._now()
    docs = list(proteins.to_documents(embeddable, vectors, source=gaborieau.SOURCE, now=now))
    n = bulk.bulk_index(
        client,
        config.PROTEINS,
        ((bulk.doc_id(gaborieau.SOURCE, pid), doc) for pid, doc in docs),
    )
    print(f"  {config.PROTEINS:20s} indexed {n} -> {bulk.count(client, config.PROTEINS):,} in index")

    # Denormalise one representative RBP vector per phage: Stage A's RBP arm is a
    # single kNN over pf-phages rather than a protein->phage join per request.
    primary = proteins.primary_rbp(embeddable)
    indexed = {
        h["_source"]["phage_id"]: h["_id"] for h in _scan_ids(client, config.PHAGES, "phage_id")
    }
    updates = [
        (
            indexed[phage_id],
            {
                "rbp_match_vector": vectors[record.protein_id].tolist(),
                "rbp_match_id": record.protein_id,
                "rbp_match_confidence": record.confidence,
            },
        )
        for phage_id, record in primary.items()
        if phage_id in indexed and record.protein_id in vectors
    ]
    applied = bulk.bulk_update(client, config.PHAGES, iter(updates))
    print(f"  attached rbp_match_vector to {applied} of {len(indexed)} phages")

    without = sorted(set(indexed) - set(primary))
    if without:
        print(f"  {len(without)} phages have no usable RBP: {without[:8]}")
    return 0


def _cosine_spread(vectors: dict) -> tuple[float, float, float]:
    """(min, max, median) pairwise cosine -- the check that centring worked."""
    import numpy as np

    matrix = np.stack(list(vectors.values()))
    similarities = matrix @ matrix.T
    upper = similarities[np.triu_indices(len(matrix), k=1)]
    return float(upper.min()), float(upper.max()), float(np.median(upper))


def _scan_ids(client, index: str, field: str):
    """Yield all documents in an index with just ``field`` in ``_source``."""
    result = client.search(
        index=index, size=10_000, _source=[field], query={"match_all": {}}
    )
    return result["hits"]["hits"]


def cmd_funnel(args: argparse.Namespace) -> int:
    """Run the full funnel for one strain (Phase 4)."""
    from phageforge.funnel import pipeline

    client = _client()
    mappings.create_indices(client, (config.PREDICTIONS,))

    run = pipeline.run_funnel(
        client,
        args.strain,
        top_n=args.top_n,
        neighbour_k=args.neighbours,
        use_rbp=not args.no_rbp,
        use_prior=not args.no_prior,
        explain=not args.no_explain,
        build_cocktail=not args.no_cocktail,
    )
    print()
    print(pipeline.format_run(run))
    return 0 if run.ok else 2


def cmd_bench(args: argparse.Namespace) -> int:
    """Run the strain-grouped benchmark (Phase 6)."""
    from phageforge.bench import harness

    client = _client()
    print(f"\nbenchmark: {args.folds}-fold strain-grouped CV")
    results = harness.run_benchmark(
        client,
        n_folds=args.folds,
        limit=args.limit,
        neighbour_k=args.neighbours,
        methods=args.methods.split(",") if args.methods else None,
    )
    print()
    print(harness.format_table(results))
    path = harness.save(results, args.out)
    print(f"\nwritten to {path}  ({results['wall_seconds']}s)")
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    client = _client()
    print()
    for index in config.LOCAL_INDICES:
        if client.indices.exists(index=index):
            print(f"  {index:20s} {bulk.count(client, index):>10,} docs")
        else:
            print(f"  {index:20s} {'(missing)':>10}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="phageforge")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("indices", help="create index mappings")
    p.add_argument("--recreate", action="store_true", help="delete and recreate (destructive)")
    p.set_defaults(func=cmd_indices)

    p = sub.add_parser("ingest", help="parse and index the Gaborieau dataset")
    p.add_argument("--recreate", action="store_true", help="delete and recreate (destructive)")
    p.set_defaults(func=cmd_ingest)

    p = sub.add_parser("features", help="compute and attach genome vectors")
    p.set_defaults(func=cmd_features)

    p = sub.add_parser("proteins", help="recover RBPs and embed them with ESM-2")
    p.add_argument(
        "--resolve-only",
        action="store_true",
        help="report the locus-tag resolution and stop before embedding",
    )
    p.add_argument("--batch-size", type=int, default=4, help="ESM-2 batch size (CPU)")
    p.add_argument(
        "--recreate",
        action="store_true",
        help="delete and recreate pf-proteins before indexing",
    )
    p.set_defaults(func=cmd_proteins)

    p = sub.add_parser("funnel", help="run the funnel for one strain")
    p.add_argument("--strain", required=True, help="strain_id, e.g. ECOR-54")
    p.add_argument("--top-n", type=int, default=config.TOP_N)
    p.add_argument("--neighbours", type=int, default=config.NEIGHBOUR_K)
    p.add_argument("--no-rbp", action="store_true", help="ablation: drop Stage A's RBP arm")
    p.add_argument("--no-prior", action="store_true", help="ablation: zero the Stage B prior")
    p.add_argument("--no-explain", action="store_true", help="skip Stage C")
    p.add_argument("--no-cocktail", action="store_true", help="skip Stage F")
    p.set_defaults(func=cmd_funnel)

    p = sub.add_parser("bench", help="run the strain-grouped benchmark")
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--limit", type=int, help="evaluate only the first N strains (smoke test)")
    p.add_argument("--neighbours", type=int, default=config.NEIGHBOUR_K)
    p.add_argument("--methods", help="comma-separated subset, e.g. funnel,generalist")
    p.add_argument("--out", type=Path, help="where to write benchmark.json")
    p.set_defaults(func=cmd_bench)

    p = sub.add_parser("info", help="show index document counts")
    p.set_defaults(func=cmd_info)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
