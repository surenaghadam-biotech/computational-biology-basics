"""Session 16 Challenge: HPV-16 Proteome Miner.

This script parses a complete GenBank genome record (NC_001526.4), extracts
all functional Coding Sequences (CDS), builds validated SeqRecord objects,
exports them into a Multi-FASTA file, and generates comparative proteome
metrics.
"""

import os
from typing import List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_cds_data(gb_path: str) -> List[SeqRecord]:
    """Parse a GenBank file and extract valid CDS features into SeqRecord objects.

    Args:
        gb_path: Path to the input GenBank file.

    Returns:
        List of SeqRecord objects representing extracted viral proteins.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If no valid CDS features with translation are found.
    """
    if not os.path.exists(gb_path):
        raise FileNotFoundError(f"GenBank file not found: {gb_path}")

    record = SeqIO.read(gb_path, "genbank")
    extracted_proteins: List[SeqRecord] = []

    for feature in record.features:
        if feature.type == "CDS":
            qualifiers = feature.qualifiers

            # Ensure the CDS feature has an explicit amino acid translation
            if "translation" not in qualifiers:
                continue

            amino_acids = qualifiers["translation"][0]

            # Extract protein ID and gene symbol safely
            protein_id = qualifiers.get(
                "protein_id", [qualifiers.get("gene", ["unknown_protein"])[0]]
            )[0]
            gene_symbol = qualifiers.get("gene", ["N/A"])[0]
            product_desc = qualifiers.get("product", ["hypothetical protein"])[
                0
            ]

            protein_record = SeqRecord(
                seq=Seq(amino_acids),
                id=protein_id,
                name=gene_symbol,
                description=f"gene={gene_symbol} | product={product_desc}",
            )
            extracted_proteins.append(protein_record)

    if not extracted_proteins:
        raise ValueError(f"No translated CDS features found in {gb_path}")

    return extracted_proteins


def export_proteome_fasta(
    protein_list: List[SeqRecord], output_path: str
) -> None:
    """Write extracted protein records to a Multi-FASTA file.

    Args:
        protein_list: List of protein SeqRecords.
        output_path: File system path for the output FASTA file.

    Raises:
        ValueError: If the protein list is empty.
    """
    if not protein_list:
        raise ValueError("Protein list is empty. Nothing to export.")

    count = SeqIO.write(protein_list, output_path, "fasta")
    print(f"[SUCCESS] Exported {count} protein sequences to: {output_path}")


def display_proteome_summary(protein_list: List[SeqRecord]) -> None:
    """Compute and display summary metrics for the extracted proteome.

    Args:
        protein_list: List of protein SeqRecords.

    Raises:
        ValueError: If the protein list is empty.
    """
    if not protein_list:
        raise ValueError("Cannot summarize an empty proteome.")

    # Using max and min with lambda for clean, pythonic record comparison
    max_protein = max(protein_list, key=lambda rec: len(rec.seq))
    min_protein = min(protein_list, key=lambda rec: len(rec.seq))

    print("\n" + "=" * 60)
    print("           HPV-16 PROTEOME ANALYTICS REPORT")
    print("=" * 60)
    print(f"Total Extracted CDS/Proteins : {len(protein_list)}")
    print(
        f"Longest Protein              : {max_protein.name} ({max_protein.id}) -> {len(max_protein.seq)} aa"
    )
    print(
        f"Shortest Protein             : {min_protein.name} ({min_protein.id}) -> {len(min_protein.seq)} aa"
    )
    print("-" * 60)
    print("Extracted Gene Catalog:")
    for idx, p in enumerate(protein_list, 1):
        print(
            f"  [{idx:02d}] {p.name:<6} | {p.id:<14} | Length: {len(p.seq):>4} aa"
        )
    print("=" * 60 + "\n")


if __name__ == "__main__":
    INPUT_GB = "hpv16_genome.gb"
    OUTPUT_FASTA = "hpv16_proteome.fasta"

    try:
        proteome = extract_cds_data(INPUT_GB)
        export_proteome_fasta(proteome, OUTPUT_FASTA)
        display_proteome_summary(proteome)
    except Exception as exc:
        print(f"[ERROR] Pipeline execution failed: {exc}")
