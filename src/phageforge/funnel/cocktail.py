"""Stage F -- cocktail design by greedy set-cover with a receptor-diversity penalty.

The only scoring the architecture permits in Python (TECHNICAL_DESIGN.md R3),
because it is a combinatorial optimisation over an already-ranked shortlist
rather than a retrieval or ranking operation.

The biology: a single receptor mutation escapes *every* phage that binds that
receptor, so four phages that all target the same LPS core are one phage wearing
four hats. The objective therefore penalises within-set receptor overlap:

    score(S) = sum_p P(infect | p, strain) - lambda * sum_{p != q} overlap(p, q)

**Receptor overlap is a proxy here, and must be labelled as one.** The designed
key is ``receptor_classes`` (ARCHITECTURE.md 5.1), derived from RBP homology and
annotation keywords -- which this archive does not ship (DATA_NOTES.md). What we
do have is the ESM-2 embedding of each phage's primary RBP, so overlap is
measured as cosine similarity between those embeddings: two phages whose fibres
are near-identical in ESM-2 space almost certainly engage the same receptor.
That is a defensible stand-in and a genuinely useful one, but it is *inferred*
similarity, not an assigned receptor class, and the UI should say so.
"""

from __future__ import annotations

import math
from collections.abc import Sequence

from elasticsearch import Elasticsearch

from phageforge import config

#: Receptor-diversity pressure. At 1.0 a second phage whose RBP is identical to
#: an already-chosen one contributes nothing; at 0 the cocktail is just the
#: top-N. Exposed so a reviewer can watch the set change as it moves
#: (TECHNICAL_DESIGN.md 8.3).
LAMBDA = 1.0

#: Cocktails larger than this stop being a practical bench protocol.
MAX_SIZE = 4


def _cosine(a: Sequence[float], b: Sequence[float]) -> float:
    dot = sum(x * y for x, y in zip(a, b, strict=True))
    na = math.sqrt(sum(x * x for x in a))
    nb = math.sqrt(sum(y * y for y in b))
    if na == 0 or nb == 0:
        return 0.0
    return dot / (na * nb)


def design(
    client: Elasticsearch,
    shortlist: Sequence[dict],
    *,
    size: int = MAX_SIZE,
    lam: float = LAMBDA,
) -> list[dict]:
    """Greedy set-cover over a finished shortlist. Returns the chosen phages.

    Each entry records the marginal gain it contributed and why it was taken, so
    the cocktail can be explained rather than just displayed.
    """
    if not shortlist:
        return []

    ids = [p["phage_id"] for p in shortlist]
    hits = client.search(
        index=config.PHAGES,
        size=len(ids),
        query={"terms": {"phage_id": ids}},
        source_includes=["phage_id", "rbp_match_vector", "rbp_match_confidence"],
    )["hits"]["hits"]
    vectors = {
        h["_source"]["phage_id"]: h["_source"].get("rbp_match_vector")
        for h in hits
    }

    # Probability proxy: the neighbour-transfer prior when there is one, else the
    # normalised Stage D score. Stated as a proxy -- these are not calibrated
    # probabilities and the API labels them accordingly.
    top_score = max((p["score"] for p in shortlist), default=1.0) or 1.0
    probability = {
        p["phage_id"]: p["prior"] if p.get("n_neighbours_tested") else p["score"] / top_score
        for p in shortlist
    }

    def overlap(a: str, b: str) -> float:
        va, vb = vectors.get(a), vectors.get(b)
        if not va or not vb:
            # No RBP for one of them: assume no known overlap, but the entry
            # below records that the diversity claim is unverified.
            return 0.0
        return max(0.0, _cosine(va, vb))

    chosen: list[dict] = []
    remaining = list(ids)
    while remaining and len(chosen) < size:
        best_id, best_gain, best_overlap = None, None, 0.0
        for candidate in remaining:
            penalty = sum(overlap(candidate, c["phage_id"]) for c in chosen)
            gain = probability.get(candidate, 0.0) - lam * penalty
            if best_gain is None or gain > best_gain:
                best_id, best_gain, best_overlap = candidate, gain, penalty
        if best_id is None:
            break
        # Stop once adding another phage actively hurts: a cocktail of near-clones
        # is worse than a smaller diverse one.
        if chosen and best_gain is not None and best_gain <= 0:
            break

        has_vector = bool(vectors.get(best_id))
        if not chosen:
            rationale = "highest-probability phage, seeds the cocktail"
        elif not has_vector:
            rationale = "no RBP embedding, receptor diversity unverified"
        elif best_overlap < 0.5:
            rationale = f"RBP distinct from the set (max overlap {best_overlap:.2f})"
        else:
            rationale = f"RBP similar to the set (overlap {best_overlap:.2f}), added on probability"

        chosen.append(
            {
                "phage_id": best_id,
                "probability_proxy": round(probability.get(best_id, 0.0), 4),
                "receptor_overlap": round(best_overlap, 4),
                "marginal_gain": round(best_gain or 0.0, 4),
                "rbp_confidence": next(
                    (h["_source"].get("rbp_match_confidence") for h in hits
                     if h["_source"]["phage_id"] == best_id),
                    None,
                ),
                "overlap_basis": "esm2_rbp_cosine" if has_vector else "none",
                "rationale": rationale,
            }
        )
        remaining.remove(best_id)
    return chosen
