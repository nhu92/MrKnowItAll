from pathlib import Path
from zipfile import ZipFile

import pytest

from sprout_refbuilder.backbone import Backbone, install_backbone_archive


def test_install_backbone_archive(tmp_path: Path) -> None:
    archive = tmp_path / "backbone.zip"
    with ZipFile(archive, "w") as bundle:
        bundle.writestr("nested/1.fasta", ">x\nAC-G\n>y\nA--G\n")
        bundle.writestr("nested/gene_summary.csv", "Gene,Count\n1,2\n")
        bundle.writestr("nested/ignored.py", "raise RuntimeError\n")

    destination = tmp_path / "installed"
    report = install_backbone_archive(archive, destination, expected_loci=1)

    assert report["locus_count"] == 1
    assert report["sequence_count"] == 2
    assert report["source_sha256"]
    assert (destination / "1.fasta").is_file()
    assert (destination / "gene_summary.csv").is_file()
    assert not (destination / "ignored.py").exists()
    assert Backbone(destination).metadata()["locus_count"] == 1

    with pytest.raises(FileExistsError):
        install_backbone_archive(archive, destination, expected_loci=1)


def test_install_rejects_duplicate_flattened_names(tmp_path: Path) -> None:
    archive = tmp_path / "duplicate.zip"
    with ZipFile(archive, "w") as bundle:
        bundle.writestr("a/1.fasta", ">x\nAC\n")
        bundle.writestr("b/1.fasta", ">y\nAC\n")

    with pytest.raises(ValueError, match="duplicate"):
        install_backbone_archive(archive, tmp_path / "installed", expected_loci=2)
