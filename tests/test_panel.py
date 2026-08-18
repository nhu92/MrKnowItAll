from pathlib import Path

from sprout_refbuilder.panel import HierarchicalPanelBuilder


class FakeKewClient:
    def __init__(self, tree: Path) -> None:
        self.tree = tree

    def sync(self) -> dict:
        return {"release": "test"}

    def fetch_species_tree(self) -> Path:
        return self.tree


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_family_panel_uses_current_kew_taxonomy(tmp_path: Path) -> None:
    backbone = tmp_path / "backbone"
    sequence = "ATG" + "ACG" * 20
    write(
        backbone / "4691.fasta",
        f">Asparagales_Alliaceae_Allium_sativum_1KP-4691\n{sequence}\n"
        f">Asparagales_Asparagaceae_Asparagus_officinalis_JGI-4691\n{sequence}\n"
        f">Rosales_Moraceae_Ficus_carica_1KP-4691\n{sequence}\n",
    )
    summary = tmp_path / "species_summary.csv"
    write(
        summary,
        "Species,Order,Family,Data Source,Gene Count\n"
        "Allium_sativum,Asparagales,Alliaceae,1KP,300\n"
        "Asparagus_officinalis,Asparagales,Asparagaceae,JGI,250\n"
        "Ficus_carica,Rosales,Moraceae,1KP,200\n",
    )
    tree = tmp_path / "tree.nwk"
    write(
        tree,
        "(Asparagales_Asparagaceae_Allium_sativum_ERR1,"
        "Asparagales_Asparagaceae_Asparagus_officinalis_ERR2,"
        "Rosales_Moraceae_Ficus_carica_ERR3);",
    )
    output = tmp_path / "panel"
    report = HierarchicalPanelBuilder(
        backbone, summary, output, representatives_per_taxon=2, prefer_mafft=False
    ).build("family", FakeKewClient(tree), parent_taxa=["Asparagales"])

    result = (output / "ref" / "4691.fasta").read_text()
    assert report["selected_groups"] == 1
    assert report["selected_specimens"] == 2
    assert result.count(">") == 2
    assert "Asparagales_Asparagaceae_Allium_sativum" in result
    assert "Ficus" not in result
