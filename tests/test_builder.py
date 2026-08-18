import json
from pathlib import Path

from sprout_refbuilder.builder import ReferenceBuilder
from sprout_refbuilder.models import Taxonomy


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_build_local_reference_bundle(tmp_path: Path) -> None:
    backbone = tmp_path / "backbone"
    sequence = "ATG" + "ACG" * 40 + "TAA"
    write(
        backbone / "4691.fasta",
        f">Brassicales_Salicaceae_Other_species_A\n{sequence}\n"
        f">Brassicales_Salicaceae_Another_species_B\n{sequence}\n",
    )
    recovery = tmp_path / "recovery.fasta"
    write(recovery, f">4691 Gene_Name:test Species:Abatia_rugosa\n{sequence}\n")
    output = tmp_path / "result"

    report = ReferenceBuilder(backbone, output, prefer_mafft=False).build_from_fasta(
        "Abatia rugosa",
        Taxonomy("Malpighiales", "Salicaceae", "Abatia", "rugosa"),
        recovery,
        "ERR1",
    )

    aligned = (output / "ref" / "4691.fasta").read_text()
    assert report.accepted_loci == ["4691"]
    assert aligned.count(">") == 3
    assert "Malpighiales_Salicaceae_Abatia_rugosa_ERR1" in aligned
    assert Path(str(output) + "_bundle.zip").exists()
    import zipfile

    with zipfile.ZipFile(Path(str(output) + "_bundle.zip")) as archive:
        assert "report.json" in archive.namelist()
    saved = json.loads((output / "report.json").read_text())
    assert saved["provenance"]["sprout_reference_argument"].endswith("ref")
