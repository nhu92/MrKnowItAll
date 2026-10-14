from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Lock
from time import sleep

import httpx
import pytest

from sprout_refbuilder.kew import KewClient
from sprout_refbuilder.models import KewRecord


def test_sync_search_and_fetch(tmp_path: Path) -> None:
    responses = {
        "/sequence_manifest.txt": (
            "Repository_name\tSequence_ID\tSequence_type\tSpecies_name\tProject_name\n"
            "INSDC\tERR1\tRead\tAbatia rugosa\tPAFTOL\n"
            "INSDC\tGCA_01.2\tGenome\tTesta species\tINSDC\n"
        ),
        "/fasta/by_recovery/": (
            '<a href="INSDC.ERR1.Abatia_rugosa.a353.fasta">x</a>'
            '<a href="INSDC.GCA_01.2.Testa_species.a353.fasta">y</a>'
        ),
        "/": '<a href="kew_tree_of_life_release_notes_4.0.txt">notes</a>',
        "/kew_tree_of_life_release_notes_4.0.txt": "Kew Tree of Life Explorer Release 4.0",
        "/fasta/by_recovery/INSDC.ERR1.Abatia_rugosa.a353.fasta": ">4691 Gene_Name:x\nACGT\n",
    }

    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, text=responses[request.url.path])

    with KewClient(tmp_path, "https://example.test/", transport=httpx.MockTransport(handler)) as client:
        metadata = client.sync()
        match = client.resolve_exact("Abatia_rugosa")[0]
        genome = client.resolve_exact("Testa species")[0]
        path = client.fetch_recovery(match)

    assert metadata["release"] == "4.0"
    assert match.recovery_file == "INSDC.ERR1.Abatia_rugosa.a353.fasta"
    assert genome.recovery_file == "INSDC.GCA_01.2.Testa_species.a353.fasta"
    assert path.read_text().startswith(">4691")


def test_offline_cache_miss_fails_without_http(tmp_path: Path) -> None:
    with (
        KewClient(tmp_path, "https://example.test/", offline=True) as client,
        pytest.raises(RuntimeError, match="Kew cache miss in offline mode"),
    ):
        client.sync()


def test_offline_recovery_miss_fails_without_http(tmp_path: Path) -> None:
    record = KewRecord("INSDC", "ERR1", "read", "Abatia rugosa", "PAFTOL", "a.fasta")
    with (
        KewClient(tmp_path, "https://example.test/", offline=True) as client,
        pytest.raises(RuntimeError, match="Kew recovery is missing"),
    ):
        client.fetch_recovery(record)


def test_concurrent_metadata_initialization_is_complete(tmp_path: Path) -> None:
    responses = {
        "/sequence_manifest.txt": (
            "Repository_name\tSequence_ID\tSequence_type\tSpecies_name\tProject_name\n"
            "INSDC\tERR1\tRead\tAbatia rugosa\tPAFTOL\n"
        ),
        "/fasta/by_recovery/": '<a href="INSDC.ERR1.Abatia_rugosa.a353.fasta">x</a>',
        "/": '<a href="kew_tree_of_life_release_notes_4.0.txt">notes</a>',
        "/kew_tree_of_life_release_notes_4.0.txt": "Release 4.0",
        "/tree/species/treeoflife.current.tree": "(Abatia_rugosa);",
    }
    calls: list[str] = []
    call_lock = Lock()

    def handler(request: httpx.Request) -> httpx.Response:
        with call_lock:
            calls.append(request.url.path)
        sleep(0.01)
        return httpx.Response(200, text=responses[request.url.path])

    def read_one(_: int) -> tuple[int, str]:
        with KewClient(
            tmp_path, "https://example.test/", transport=httpx.MockTransport(handler)
        ) as client:
            return len(client.records()), client.fetch_species_tree().read_text()

    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(read_one, range(4)))
    assert results == [(1, "(Abatia_rugosa);")] * 4
    assert calls.count("/sequence_manifest.txt") == 1
    assert calls.count("/tree/species/treeoflife.current.tree") == 1
