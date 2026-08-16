"""Fetch and verify the Gaborieau et al. (2024) Zenodo archive.

Record 13831957 accompanies *Prediction of strain-level phage-host interactions
across the Escherichia genus* (Nature Microbiology, 2024) and mirrors the
``mdmparis/coli_phage_interactions_2023`` repository at publication time.

The archive is the only Tier-1 source with **negative** interactions, which is
what makes ranking possible at all -- so integrity is checked, not assumed.
"""

from __future__ import annotations

import hashlib
import tarfile
from dataclasses import dataclass
from pathlib import Path

import requests

from phageforge.config import RAW_DIR

RECORD_ID = "13831957"
ARCHIVE_NAME = "coli_phage_interactions_2023.tar.gz"
ARCHIVE_URL = f"https://zenodo.org/records/{RECORD_ID}/files/{ARCHIVE_NAME}?download=1"
EXPECTED_MD5 = "380d288b8cd7629660603447b9856fc9"
EXPECTED_BYTES = 146_786_594

ARCHIVE_PATH = RAW_DIR / ARCHIVE_NAME
EXTRACT_DIR = RAW_DIR / "gaborieau"


@dataclass
class FetchResult:
    path: Path
    md5: str
    size: int
    verified: bool
    downloaded: bool


def md5sum(path: Path, chunk: int = 1 << 20) -> str:
    digest = hashlib.md5()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            digest.update(block)
    return digest.hexdigest()


def download(force: bool = False) -> FetchResult:
    """Download the archive if absent, then verify size and checksum."""
    downloaded = False

    if force or not ARCHIVE_PATH.exists() or ARCHIVE_PATH.stat().st_size != EXPECTED_BYTES:
        ARCHIVE_PATH.parent.mkdir(parents=True, exist_ok=True)
        with requests.get(ARCHIVE_URL, stream=True, timeout=120) as response:
            response.raise_for_status()
            with ARCHIVE_PATH.open("wb") as fh:
                for chunk in response.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
        downloaded = True

    size = ARCHIVE_PATH.stat().st_size
    digest = md5sum(ARCHIVE_PATH)
    return FetchResult(
        path=ARCHIVE_PATH,
        md5=digest,
        size=size,
        verified=(digest == EXPECTED_MD5 and size == EXPECTED_BYTES),
        downloaded=downloaded,
    )


def extract(dest: Path = EXTRACT_DIR, *, force: bool = False) -> Path:
    """Extract the archive. Refuses to run on an unverified download."""
    if dest.exists() and not force and any(dest.iterdir()):
        return dest

    result = download()
    if result.size != EXPECTED_BYTES:
        raise RuntimeError(
            f"Archive is {result.size} bytes, expected {EXPECTED_BYTES}. "
            "Download is incomplete -- re-run with force=True."
        )

    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(result.path, "r:gz") as tar:
        # filter='data' rejects absolute paths, parent traversal and special files.
        tar.extractall(dest, filter="data")
    return dest


def manifest(root: Path = EXTRACT_DIR, limit: int | None = None) -> list[tuple[str, int]]:
    """List extracted files as ``(relative_path, size_bytes)``, largest first."""
    if not root.exists():
        raise FileNotFoundError(f"{root} does not exist -- run extract() first")
    entries = [
        (str(p.relative_to(root)), p.stat().st_size)
        for p in root.rglob("*")
        if p.is_file()
    ]
    entries.sort(key=lambda kv: kv[1], reverse=True)
    return entries[:limit] if limit else entries
