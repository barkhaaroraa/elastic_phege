"""Batched bulk indexing with deterministic IDs and loud failures.

Two properties matter here:

* **Idempotency.** Document ``_id`` is ``sha1(source + native_id)``, so re-running
  ingest upserts rather than duplicating. ``make ingest`` twice must leave the
  same doc counts.
* **Errors surface.** A partially-failed bulk request returns HTTP 200 with
  per-item errors. Silently ignoring those is how you end up benchmarking against
  half a dataset, so we collect and raise.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable, Iterator
from typing import Any

from elasticsearch import Elasticsearch

BATCH_SIZE = 5_000


def doc_id(source: str, native_id: str) -> str:
    """Deterministic document ID, so re-ingest upserts instead of duplicating."""
    return hashlib.sha1(f"{source}::{native_id}".encode()).hexdigest()


def _batched(it: Iterable[dict], size: int) -> Iterator[list[dict]]:
    batch: list[dict] = []
    for item in it:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch


class BulkIndexError(RuntimeError):
    """Raised when any document in a bulk request fails to index."""

    def __init__(self, failures: list[dict]) -> None:
        self.failures = failures
        sample = failures[:3]
        super().__init__(
            f"{len(failures)} document(s) failed to index. First failures: {sample}"
        )


def bulk_index(
    client: Elasticsearch,
    index: str,
    docs: Iterable[tuple[str, dict[str, Any]]],
    *,
    batch_size: int = BATCH_SIZE,
    refresh_at_end: bool = True,
) -> int:
    """Index ``(doc_id, source_doc)`` pairs. Returns the number indexed.

    Refresh is disabled during the load and performed once at the end -- refreshing
    per batch is the usual cause of a slow bulk ingest.
    """
    total = 0
    failures: list[dict] = []

    for batch in _batched(({"_id": i, "_doc": d} for i, d in docs), batch_size):
        operations: list[dict] = []
        for item in batch:
            operations.append({"index": {"_index": index, "_id": item["_id"]}})
            operations.append(item["_doc"])

        response = client.bulk(operations=operations, refresh=False)

        if response.get("errors"):
            for entry in response["items"]:
                result = entry.get("index") or entry.get("create") or {}
                if result.get("error"):
                    failures.append(
                        {"id": result.get("_id"), "error": result["error"]}
                    )
        total += len(batch)

    if failures:
        raise BulkIndexError(failures)

    if refresh_at_end:
        client.indices.refresh(index=index)
    return total


def bulk_update(
    client: Elasticsearch,
    index: str,
    updates: Iterable[tuple[str, dict[str, Any]]],
    *,
    batch_size: int = BATCH_SIZE,
    refresh_at_end: bool = True,
) -> int:
    """Apply partial ``(doc_id, partial_doc)`` updates. Returns the count applied.

    Used to attach vectors to documents that were already ingested, rather than
    rebuilding and re-sending the full source.
    """
    total = 0
    failures: list[dict] = []

    for batch in _batched(({"_id": i, "_doc": d} for i, d in updates), batch_size):
        operations: list[dict] = []
        for item in batch:
            operations.append({"update": {"_index": index, "_id": item["_id"]}})
            operations.append({"doc": item["_doc"]})

        response = client.bulk(operations=operations, refresh=False)
        if response.get("errors"):
            for entry in response["items"]:
                result = entry.get("update") or {}
                if result.get("error"):
                    failures.append({"id": result.get("_id"), "error": result["error"]})
        total += len(batch)

    if failures:
        raise BulkIndexError(failures)

    if refresh_at_end:
        client.indices.refresh(index=index)
    return total


def count(client: Elasticsearch, index: str, query: dict | None = None) -> int:
    """Document count, optionally filtered."""
    if query is None:
        return int(client.count(index=index)["count"])
    return int(client.count(index=index, query=query)["count"])
