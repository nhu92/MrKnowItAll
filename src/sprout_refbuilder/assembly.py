from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


def assemble_with_hybpiper(
    read1: Path,
    read2: Path,
    target_file: Path,
    taxon_slug: str,
    work_dir: Path,
    threads: int = 8,
) -> Path:
    """Run HybPiper and return its retrieved nucleotide FASTA.

    This is deliberately isolated from Kew assembled-recovery mode. The target file is only
    used by HybPiper for read recovery; it is never emitted as a SPrOUT reference alignment.
    """
    if not shutil.which("hybpiper"):
        raise RuntimeError("hybpiper is required for local FASTQ mode but was not found on PATH")
    work_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            "hybpiper",
            "assemble",
            "-t_dna",
            str(target_file.resolve()),
            "-r",
            str(read1.resolve()),
            str(read2.resolve()),
            "--prefix",
            taxon_slug,
            "--bwa",
            "--cpu",
            str(threads),
        ],
        cwd=work_dir,
        check=True,
    )
    name_list = work_dir / "namelist.txt"
    name_list.write_text(taxon_slug + "\n", encoding="utf-8")
    subprocess.run(
        [
            "hybpiper",
            "retrieve_sequences",
            "dna",
            str(target_file.resolve()),
            str(name_list),
        ],
        cwd=work_dir,
        check=True,
    )
    expected = work_dir / f"{taxon_slug}.fasta"
    if not expected.exists():
        matches = sorted(work_dir.glob("*.FNA")) + sorted(work_dir.glob("*.fasta"))
        if not matches:
            raise RuntimeError("HybPiper completed but no retrieved nucleotide FASTA was found")
        expected = matches[0]
    return expected

