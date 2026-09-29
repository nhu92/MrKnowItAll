import csv
from pathlib import Path

from sprout_refbuilder.candidates import select_mixture_candidates


def _write_matrix(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        ",mix_query_a,mix_query_b,"
        "Rosales_FamA_GenusA_one_X,Rosales_FamB_GenusB_two_Y\n"
        "mix_query_a,0,1,0.1,0.9\n"
        "mix_query_b,1,0,0.9,0.1\n"
        "Rosales_FamA_GenusA_one_X,0.1,0.9,0,1\n"
        "Rosales_FamB_GenusB_two_Y,0.9,0.1,1,0\n",
        encoding="utf-8",
    )


def test_locus_votes_rescue_low_aggregate_mixture_member(tmp_path: Path) -> None:
    matrix_dir = tmp_path / "matrix"
    _write_matrix(matrix_dir / "gene1.cleaned.csv")
    _write_matrix(matrix_dir / "gene2.cleaned.csv")
    predictions = tmp_path / "predictions.csv"
    predictions.write_text(
        "row_name,sum_of_total_value,z_score\nFamA,100,1.0\nFamB,1,-1.0\n",
        encoding="utf-8",
    )
    output = tmp_path / "candidates.txt"
    report = tmp_path / "evidence.csv"

    result = select_mixture_candidates(
        matrix_dir,
        predictions,
        output,
        report,
        level="family",
        project="mix",
        minimum_locus_votes=2,
    )

    assert result["selected"] == ["FamA", "FamB"]
    assert output.read_text(encoding="utf-8").splitlines() == ["FamA", "FamB"]
    with report.open(encoding="utf-8", newline="") as handle:
        rows = {row["taxon"]: row for row in csv.DictReader(handle)}
    assert rows["FamB"]["selected_by"] == "locus_votes"
    assert rows["FamB"]["locus_votes"] == "2"
