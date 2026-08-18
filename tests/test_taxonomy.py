from pathlib import Path

from sprout_refbuilder.taxonomy import taxonomy_from_kew_tree


def test_taxonomy_from_tree(tmp_path: Path) -> None:
    tree = tmp_path / "tree.nwk"
    tree.write_text(
        "(Malpighiales_Salicaceae_Abatia_rugosa_ERR1:0.1,"
        "Rosales_Rosaceae_Rosa_rubiginosa_ERR2:0.2);",
        encoding="utf-8",
    )
    taxonomy = taxonomy_from_kew_tree(tree, "Abatia rugosa", "ERR1")
    assert taxonomy.sprout_prefix == "Malpighiales_Salicaceae_Abatia_rugosa"
