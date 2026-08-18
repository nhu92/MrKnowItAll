from pathlib import Path

import httpx

from sprout_refbuilder.kew import KewClient


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
