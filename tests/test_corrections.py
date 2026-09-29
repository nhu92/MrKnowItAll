import csv
from pathlib import Path

from sprout_refbuilder.corrections import (
    TaxonomyCorrection,
    apply_taxonomy_corrections,
    data_source_from_identifier,
    load_taxonomy_corrections,
)


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def test_apply_replacements_exclusions_and_regenerate_summaries(tmp_path: Path) -> None:
    directory = tmp_path / "backbone"
    write(
        directory / "4691.fasta",
        ">Oldorder_Oldfamily_Genus_species_1KP-4691 source\nAC-G\n"
        ">Oldorder_Badfamily_Bad_species_Timilsena_et_al-4691 source\nA--G\n"
        ">Otherorder_Otherfamily_Other_species_JGI-4691 source\nACGG\n",
    )
    report = apply_taxonomy_corrections(
        directory,
        [
            TaxonomyCorrection(
                "replace",
                "Oldorder_Oldfamily_Genus_species",
                "Neworder_Newfamily_Newgenus_newspecies",
                "",
            ),
            TaxonomyCorrection(
                "exclude",
                "Oldorder_Badfamily_Bad_species_Timilsena_et_al",
                "",
                "mislabelled",
            ),
        ],
    )

    text = (directory / "4691.fasta").read_text(encoding="utf-8")
    assert "Neworder_Newfamily_Newgenus_newspecies_1KP-4691" in text
    assert "Bad_species" not in text
    assert report["replaced_locus_sequences"] == 1
    assert report["excluded_locus_sequences"] == 1
    with (directory / "species_summary.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert any(row["Species"] == "Newgenus_newspecies" for row in rows)
    assert not any(row["Species"] == "Bad_species" for row in rows)


def test_data_source_normalizes_timilsena() -> None:
    assert (
        data_source_from_identifier(
            "Liliales_Ripogonaceae_Ripogonum_brevifolium_Timilsena_et_al-4691"
        )
        == "Timilsena"
    )


def test_builtin_correction_inventory() -> None:
    rows, digest = load_taxonomy_corrections()
    assert len([row for row in rows if row.action == "replace"]) == 48
    assert len([row for row in rows if row.action == "exclude"]) == 9
    assert len(digest) == 64
