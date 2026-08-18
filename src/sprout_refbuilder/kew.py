from __future__ import annotations

import csv
import json
import re
from difflib import SequenceMatcher
from html.parser import HTMLParser
from pathlib import Path
from typing import Self
from urllib.parse import quote, urljoin

import httpx

from .models import KewRecord

DEFAULT_BASE_URL = "https://sftp.kew.org/pub/treeoflife/current_release/"


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
    ) -> None:
        self.cache_dir = cache_dir
        self.base_url = base_url.rstrip("/") + "/"
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
        if not force and all(
            p.exists()
            for p in (self.manifest_path, self.recovery_index_path, self.release_metadata_path)
        ):
            return json.loads(self.release_metadata_path.read_text(encoding="utf-8"))

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
        self.manifest_path.write_text(manifest, encoding="utf-8")
        self.recovery_index_path.write_text(json.dumps(files, indent=2), encoding="utf-8")
        metadata = {
            "release": release,
            "base_url": self.base_url,
            "recovery_count": len(files),
            "release_notes": notes_name,
        }
        self.release_metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
        return metadata

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
        destination.parent.mkdir(parents=True, exist_ok=True)
        response = self.http.get(urljoin(self.base_url, "fasta/by_recovery/" + quote(record.recovery_file)))
        response.raise_for_status()
        destination.write_bytes(response.content)
        return destination

    def fetch_species_tree(self) -> Path:
        destination = self.cache_dir / "treeoflife.current.tree"
        if not destination.exists():
            destination.write_bytes(self._get_bytes("tree/species/treeoflife.current.tree"))
        return destination

    def _get_text(self, relative: str) -> str:
        return self._get_bytes(relative).decode("utf-8", errors="replace")

    def _get_bytes(self, relative: str) -> bytes:
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
