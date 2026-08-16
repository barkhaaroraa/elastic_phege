"""Central configuration: paths, index names, and the Elasticsearch client factory."""

from __future__ import annotations

import os
from pathlib import Path

from elasticsearch import Elasticsearch

# ---------------------------------------------------------------- paths

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("PF_DATA_DIR", REPO_ROOT / "data"))
RAW_DIR = DATA_DIR / "raw"
DERIVED_DIR = DATA_DIR / "derived"
CACHE_DIR = DERIVED_DIR / "cache"

for _d in (RAW_DIR, DERIVED_DIR, CACHE_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- elasticsearch

ES_URL = os.environ.get("PF_ES_URL", "http://localhost:9200")
ES_API_KEY = os.environ.get("PF_ES_API_KEY")  # set when we move to Elastic Cloud

PHAGES = "pf-phages"
BACTERIA = "pf-bacteria"
INTERACTIONS = "pf-interactions"
PROTEINS = "pf-proteins"
LITERATURE = "pf-literature"
PREDICTIONS = "pf-predictions"

#: Indices created for the local backbone. ``LITERATURE`` needs an ELSER inference
#: endpoint, so it is only created once we are on Elastic Cloud.
LOCAL_INDICES = (PHAGES, BACTERIA, INTERACTIONS, PROTEINS, PREDICTIONS)

# ---------------------------------------------------------------- vectors

#: MinHash/k-mer genome sketch width. Low-dimensional, so ``int8_hnsw`` is the
#: right quantization here -- BBQ needs higher dimensionality to pay off.
SKETCH_DIMS = 256

#: ESM-2 ``esm2_t12_35M_UR50D`` embedding width. CPU-feasible on 8 cores; the
#: 650M model (1280-dim) is the later upgrade if this proves too weak.
ESM_MODEL = os.environ.get("PF_ESM_MODEL", "esm2_t12_35M_UR50D")
ESM_DIMS = int(os.environ.get("PF_ESM_DIMS", "480"))

# ---------------------------------------------------------------- funnel

#: Stage A candidate-set size. The LTR rescorer's ``window_size`` must equal this
#: so every survivor is scored (see ARCHITECTURE.md 7.4).
CANDIDATE_POOL = 500
#: Neighbour strains pulled in Stage B for susceptibility transfer.
NEIGHBOUR_K = 25
#: Final shortlist length.
TOP_N = 10


def get_client(url: str | None = None) -> Elasticsearch:
    """Build an Elasticsearch client.

    Uses an API key when ``PF_ES_API_KEY`` is set (Elastic Cloud), otherwise an
    unauthenticated local connection.
    """
    kwargs: dict = {
        "request_timeout": 60,
        "retry_on_timeout": True,
        "max_retries": 3,
    }
    if ES_API_KEY:
        kwargs["api_key"] = ES_API_KEY
    return Elasticsearch(url or ES_URL, **kwargs)


def require_healthy(client: Elasticsearch) -> dict:
    """Fail loudly if the cluster is unreachable or red."""
    health = client.cluster.health()
    if health["status"] == "red":
        raise RuntimeError(f"Elasticsearch cluster is red: {health}")
    return dict(health)
