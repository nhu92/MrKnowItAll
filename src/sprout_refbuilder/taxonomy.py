from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .models import Taxonomy


@dataclass(frozen=True)
class KewTreeTaxon:
    taxonomy: Taxonomy
    sequence_id: str

    @property
    def scientific_name(self) -> str:
        return f"{self.taxonomy.genus} {self.taxonomy.species}"


def taxonomy_from_kew_tree(tree_path: Path, scientific_name: str, sequence_id: str = "") -> Taxonomy:
    """Extract the documented Order_Family_Genus_species prefix from a Kew Newick leaf."""
    genus, epithet = _binomial(scientific_name)
    text = tree_path.read_text(encoding="utf-8")
    marker = f"_{genus}_{epithet}_"
    starts = [m.start() for m in re.finditer(re.escape(marker), text, flags=re.IGNORECASE)]
    if not starts:
        raise LookupError(f"{scientific_name!r} was not found in the Kew species tree")
    selected = starts[0]
    if sequence_id:
        for start in starts:
            tail = text[start : start + len(marker) + len(sequence_id) + 64]
            if sequence_id.casefold() in tail.casefold():
                selected = start
                break
    prefix = text[max(0, selected - 200) : selected].split("(")[-1].split(",")[-1]
    prefix = prefix.split(")")[-1]
    prefix = re.sub(r"^[^A-Za-z]+", "", prefix)
    tokens = prefix.split("_")
    if len(tokens) < 2:
        raise LookupError(f"Could not parse order/family for {scientific_name!r}")
    order, family = tokens[-2], tokens[-1]
    return Taxonomy(order=order, family=family, genus=genus, species=epithet)


def taxonomy_from_header(description: str) -> Taxonomy | None:
    tokens = description.split()[0].split("_")
    if len(tokens) < 4:
        return None
    return Taxonomy(tokens[0], tokens[1], tokens[2], tokens[3])


def read_kew_tree_taxa(tree_path: Path) -> list[KewTreeTaxon]:
    """Parse Kew's documented five-field leaf labels from a Newick species tree."""
    text = tree_path.read_text(encoding="utf-8")
    labels = re.findall(r"(?:^|[(,])([^(),;]+?)(?=[,)])", text)
    taxa: list[KewTreeTaxon] = []
    for label in labels:
        fields = label.split("_", 4)
        if len(fields) != 5:
            continue
        order, family, genus, species, sequence_id = fields
        taxa.append(KewTreeTaxon(Taxonomy(order, family, genus, species), sequence_id))
    return taxa


def kew_taxonomy_index(tree_path: Path) -> dict[str, list[KewTreeTaxon]]:
    index: dict[str, list[KewTreeTaxon]] = {}
    for taxon in read_kew_tree_taxa(tree_path):
        index.setdefault(taxon.scientific_name.casefold(), []).append(taxon)
    return index


def _binomial(name: str) -> tuple[str, str]:
    tokens = name.replace("_", " ").split()
    if len(tokens) < 2:
        raise ValueError("A binomial scientific name (Genus species) is required")
    return tokens[0], tokens[1]
