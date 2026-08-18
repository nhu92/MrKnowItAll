from sprout_refbuilder.fasta import FastaRecord, format_fasta, parse_fasta_text


def test_fasta_round_trip() -> None:
    records = parse_fasta_text(">one description\nACGT.N-\n>two\nTTAA\n")
    assert records == [
        FastaRecord("one", "one description", "ACGT-N-"),
        FastaRecord("two", "two", "TTAA"),
    ]
    assert parse_fasta_text(format_fasta(records)) == records

