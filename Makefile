.PHONY: help venv up down destroy status indices ingest features proteins funnel api bench test lint clean info

PY := .venv/bin/python
UV := $(shell command -v uv 2>/dev/null || echo "$$HOME/.local/bin/uv")

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

venv: ## Create the Python 3.12 venv and install the project
	$(UV) venv --python 3.12 .venv
	$(UV) pip install --python .venv/bin/python -e ".[dev]"

up: ## Start local Elasticsearch (podman) and wait for health
	@bash infra/elasticsearch.sh up

down: ## Stop Elasticsearch
	@bash infra/elasticsearch.sh down

destroy: ## Remove the ES container and its data volume
	@bash infra/elasticsearch.sh destroy

status: ## Show cluster health
	@bash infra/elasticsearch.sh status

indices: ## Create the index mappings
	$(PY) -m phageforge.cli indices

ingest: ## Fetch + parse the Gaborieau dataset into Elasticsearch
	$(PY) -m phageforge.cli ingest

features: ## Compute genome sketches and phylogenetic embeddings
	$(PY) -m phageforge.cli features

proteins: ## Recover RBPs (pyrodigal) and embed them with ESM-2
	$(PY) -m phageforge.cli proteins

funnel: ## Run the funnel for one strain (STRAIN=<id>)
	$(PY) -m phageforge.cli funnel --strain $(STRAIN)

api: ## Run the FastAPI service
	.venv/bin/uvicorn phageforge.api.main:app --reload --port 8000

bench: ## Run the benchmark harness
	$(PY) -m phageforge.cli bench

info: ## Show document counts per index
	$(PY) -m phageforge.cli info

test: ## Run unit tests
	.venv/bin/pytest -q

lint: ## Lint
	.venv/bin/ruff check src tests

clean: ## Remove derived artifacts (keeps raw downloads)
	rm -rf data/derived
