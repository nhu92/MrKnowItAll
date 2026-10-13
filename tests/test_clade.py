from pathlib import Path

import pytest

from sprout_refbuilder.clade import (
    BalancedKewPanelBuilder,
    _alignment_backbone,
    _load_locus,
    read_taxon_candidates,
    taxon_name,
)
from sprout_refbuilder.models import KewRecord, Taxonomy


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


class FakeKewClient:
    base_url = "https://example.test/"

    def __init__(self, root: Path) -> None:
        self.root = root
        self.tree = root / "tree.nwk"
        write(
            self.tree,
            "(Rosales_FamA_GenusA_speciesa_A1,Rosales_FamB_GenusB_speciesb_B1);\n",
        )
        self._records = [
            KewRecord("INSDC", "A1", "Read", "GenusA speciesa", "PAFTOL", "A.fasta"),
            KewRecord("INSDC", "B1", "Read", "GenusB speciesb", "PAFTOL", "B.fasta"),
        ]
        sequence = "ATG" + "ACG" * 40 + "TAA"
        write(root / "A.fasta", f">4691 Gene_Name:test\n{sequence}\n")
        write(root / "B.fasta", f">4691 Gene_Name:test\n{sequence}\n")

    def sync(self) -> dict:
        return {"release": "test"}

    def fetch_species_tree(self) -> Path:
        return self.tree

    def records(self) -> list[KewRecord]:
        return self._records

    def fetch_recovery(self, record: KewRecord) -> Path:
        assert record.recovery_file
        return self.root / record.recovery_file


def test_read_taxon_candidates_from_txt_and_csv(tmp_path: Path) -> None:
    text = tmp_path / "candidates.txt"
    write(text, "Rosales\nRosales\n# note\n")
    assert read_taxon_candidates(text) == ["Rosales"]

    csv_path = tmp_path / "predictions.csv"
    write(csv_path, "row_name,z_score\nRosales,7.3\nFabales,0.5\n")
    assert read_taxon_candidates(csv_path, minimum_z=1.0) == ["Rosales"]


def test_species_taxon_name_is_binomial() -> None:
    taxonomy = Taxonomy("Rosales", "Moraceae", "Artocarpus", "heterophyllus")
    assert taxon_name(taxonomy, "species") == "Artocarpus_heterophyllus"


def test_genus_alignment_keeps_focal_family_and_outgroups(tmp_path: Path) -> None:
    locus_path = tmp_path / "4691.fasta"
    names = [
        "Brassicales_Brassicaceae_Brassica_oleracea_1",
        "Brassicales_Brassicaceae_Brassica_rapa_2",
        "Brassicales_Capparaceae_Capparis_spinosa_3",
        "Brassicales_Capparaceae_Capparis_sepiaria_4",
        "Rosales_Rosaceae_Rosa_canina_5",
        "Rosales_Moraceae_Ficus_carica_6",
        "Malvales_Malvaceae_Hibiscus_rosa_7",
    ]
    write(locus_path, "".join(f">{name}\nATGACG\n" for name in names))
    locus = _load_locus(locus_path)
    subset = _alignment_backbone(
        locus, "genus", {"brassicales"}, {"brassicaceae"}, {"brassica"}
    )
    ids = {record.id for record in subset.records}
    assert len(subset.records) == 5
    assert names[0] in ids and names[1] in ids
    assert names[2] in ids and names[3] not in ids
    assert names[4] in ids and names[5] not in ids
    assert names[6] in ids
    assert subset.alignment_length == locus.alignment_length


@pytest.mark.parametrize("workers", [1, 2])
def test_build_balanced_family_panel(tmp_path: Path, workers: int) -> None:
    backbone = tmp_path / "backbone"
    sequence = "ATG" + "ACG" * 40 + "TAA"
    write(
        backbone / "4691.fasta",
        f">Rosales_FamA_Other_one_X\n{sequence}\n"
        f">Rosales_FamB_Another_two_Y\n{sequence}\n",
    )
    output = tmp_path / "panel"
    client = FakeKewClient(tmp_path / "kew")
    report = BalancedKewPanelBuilder(
        backbone,
        output,
        representatives_per_taxon=1,
        minimum_recovered_loci=1,
        prefer_mafft=False,
        download_workers=1,
        workers=workers,
    ).build("family", ["Rosales"], client)  # type: ignore[arg-type]

    alignment = (output / "ref" / "4691.fasta").read_text(encoding="utf-8")
    assert report["selected_groups"] == 2
    assert report["loci"] == 1
    assert report["workers"] == workers
    assert report["selection"][0]["scientific_name"].startswith("Genus")
    assert " " in report["selection"][0]["species"]
    assert alignment.count(">") == 2
    assert "Rosales_FamA_GenusA_speciesa_KEW_A1" in alignment
    assert "Rosales_FamB_GenusB_speciesb_KEW_B1" in alignment
    assert Path(str(output) + "_bundle.zip").exists()
