from __future__ import annotations

import csv
import json
import os
import re
import tempfile
import threading
from contextlib import contextmanager
from difflib import SequenceMatcher
from html.parser import HTMLParser
from pathlib import Path
from typing import Self
from urllib.parse import quote, urljoin

import httpx

from .models import KewRecord

DEFAULT_BASE_URL = "https://sftp.kew.org/pub/treeoflife/current_release/"
_LOCAL_CACHE_LOCK = threading.RLock()


@contextmanager
def _cache_metadata_lock(cache_dir: Path):
    """Coordinate Kew metadata initialization across threads and HPC processes."""
    with _LOCAL_CACHE_LOCK:
        with (cache_dir / ".metadata.lock").open("a+b") as handle:
            if os.name == "posix":
                import fcntl

                fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            try:
                yield
            finally:
                if os.name == "posix":
                    fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def _atomic_write(destination: Path, content: bytes) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp", delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(content)
    try:
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


class _LinkParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag != "a":
            return
        href = dict(attrs).get("href")
        if href:
            self.links.append(href)


class KewClient:
    """Client for Kew's documented public release tree, not an unstable private API."""

    def __init__(
        self,
        cache_dir: Path,
        base_url: str = DEFAULT_BASE_URL,
        timeout: float = 120.0,
        transport: httpx.BaseTransport | None = None,
        offline: bool = False,
    ) -> None:
        self.cache_dir = cache_dir
        self.base_url = base_url.rstrip("/") + "/"
        self.offline = offline
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.http = httpx.Client(
            timeout=timeout,
            follow_redirects=True,
            transport=transport,
            headers={"User-Agent": "sprout-refbuilder/0.1 (+https://github.com)"},
        )

    def close(self) -> None:
        self.http.close()

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    @property
    def manifest_path(self) -> Path:
        return self.cache_dir / "sequence_manifest.tsv"

    @property
    def recovery_index_path(self) -> Path:
        return self.cache_dir / "by_recovery_index.json"

    @property
    def release_metadata_path(self) -> Path:
        return self.cache_dir / "release.json"

    def sync(self, force: bool = False) -> dict:
        with _cache_metadata_lock(self.cache_dir):
            if not force:
                cached = self._complete_metadata()
                if cached is not None:
                    return cached

            manifest = self._get_text("sequence_manifest.txt")
            recovery_html = self._get_text("fasta/by_recovery/")
            root_html = self._get_text("")
            notes_name = self._find_release_notes(root_html)
            notes = self._get_text(notes_name)
            release_match = re.search(r"Release\s+(\d+(?:\.\d+)*)", notes, re.IGNORECASE)
            release = release_match.group(1) if release_match else "unknown"

            parser = _LinkParser()
            parser.feed(recovery_html)
            files = sorted(
                href for href in parser.links if href.endswith(".a353.fasta") and "/" not in href
            )
            _atomic_write(self.manifest_path, manifest.encode("utf-8"))
            _atomic_write(self.recovery_index_path, json.dumps(files, indent=2).encode("utf-8"))
            metadata = {
                "release": release,
                "base_url": self.base_url,
                "recovery_count": len(files),
                "release_notes": notes_name,
            }
            _atomic_write(
                self.release_metadata_path, json.dumps(metadata, indent=2).encode("utf-8")
            )
            return metadata

    def _complete_metadata(self) -> dict | None:
        try:
            if not all(
                path.is_file() and path.stat().st_size > 0
                for path in (self.manifest_path, self.recovery_index_path, self.release_metadata_path)
            ):
                return None
            metadata = json.loads(self.release_metadata_path.read_text(encoding="utf-8"))
            files = json.loads(self.recovery_index_path.read_text(encoding="utf-8"))
            if not isinstance(metadata, dict) or not isinstance(files, list):
                return None
            required = {"Repository_name", "Sequence_ID", "Sequence_type", "Species_name", "Project_name"}
            with self.manifest_path.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if not required.issubset(reader.fieldnames or []):
                    return None
                for row in reader:
                    if any(row.get(field) is None for field in required):
                        return None
            return metadata
        except (OSError, ValueError, TypeError, KeyError, csv.Error):
            return None

    def records(self) -> list[KewRecord]:
        self.sync()
        files = json.loads(self.recovery_index_path.read_text(encoding="utf-8"))
        # Do not split filenames on dots: genome accessions commonly contain a version suffix
        # such as GCA_043235775.1. Construct the documented name from manifest fields instead.
        file_lookup = {filename.casefold(): filename for filename in files}

        rows: list[KewRecord] = []
        with self.manifest_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                species = row["Species_name"].strip()
                repository = row["Repository_name"].strip()
                sequence_id = row["Sequence_ID"].strip()
                expected = f"{repository}.{sequence_id}.{species.replace(' ', '_')}.a353.fasta"
                filename = file_lookup.get(expected.casefold())
                rows.append(
                    KewRecord(
                        repository=repository,
                        sequence_id=sequence_id,
                        sequence_type=row["Sequence_type"].strip(),
                        species=species,
                        project=row["Project_name"].strip(),
                        recovery_file=filename,
                    )
                )
        return rows

    def search(self, query: str, limit: int = 20) -> list[KewRecord]:
        needle = _normalize_taxon(query)
        candidates = []
        for record in self.records():
            name = _normalize_taxon(record.species)
            if needle in name or name in needle:
                rank = (name != needle, abs(len(name) - len(needle)), name)
                candidates.append((rank, record))
        if not candidates:
            candidates = [
                ((True, 1.0 - SequenceMatcher(None, needle, _normalize_taxon(r.species)).ratio(), r.species), r)
                for r in self.records()
            ]
        candidates.sort(key=lambda item: item[0])
        return [record for _, record in candidates[:limit]]

    def resolve_exact(self, taxon: str) -> list[KewRecord]:
        needle = _normalize_taxon(taxon)
        return [r for r in self.records() if _normalize_taxon(r.species) == needle and r.recovery_file]

    def fetch_recovery(self, record: KewRecord) -> Path:
        if not record.recovery_file:
            raise ValueError(f"No assembled recovery is indexed for {record.species} ({record.sequence_id})")
        destination = self.cache_dir / "recoveries" / record.recovery_file
        if destination.exists() and destination.stat().st_size > 0:
            return destination
        if self.offline:
            raise RuntimeError(f"Kew recovery is missing from the offline cache: {destination}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        response = self.http.get(urljoin(self.base_url, "fasta/by_recovery/" + quote(record.recovery_file)))
        response.raise_for_status()
        _atomic_write(destination, response.content)
        return destination

    def fetch_species_tree(self) -> Path:
        destination = self.cache_dir / "treeoflife.current.tree"
        if destination.is_file() and destination.stat().st_size > 0:
            return destination
        with _cache_metadata_lock(self.cache_dir):
            if not destination.is_file() or destination.stat().st_size == 0:
                _atomic_write(destination, self._get_bytes("tree/species/treeoflife.current.tree"))
        return destination

    def _get_text(self, relative: str) -> str:
        return self._get_bytes(relative).decode("utf-8", errors="replace")

    def _get_bytes(self, relative: str) -> bytes:
        if self.offline:
            raise RuntimeError(
                f"Kew cache miss in offline mode: {urljoin(self.base_url, relative)}"
            )
        response = self.http.get(urljoin(self.base_url, relative))
        response.raise_for_status()
        return response.content

    @staticmethod
    def _find_release_notes(root_html: str) -> str:
        match = re.search(r'href="(kew_tree_of_life_release_notes_[^"]+\.txt)"', root_html)
        if not match:
            raise RuntimeError("Kew release notes were not found in the current release index")
        return match.group(1)


def _normalize_taxon(value: str) -> str:
    return " ".join(value.replace("_", " ").split()).casefold()
