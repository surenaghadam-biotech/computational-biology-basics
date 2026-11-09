"""Session 17: Genomic Feature Extraction and Viral Integration Analysis.

This module provides tools to inspect genomic features, resolve compound
spliced coordinates, extract raw nucleotide sequences directly from FeatureLocation
objects using Bio.SeqFeature, and evaluate oncogenic gene disruption events
modeled after cervical carcinogenesis mechanisms (ZJU/PubMed).
"""

import os
from typing import Dict, List, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


def load_genomic_record(gb_path: str) -> SeqRecord:
    """Load and validate a single GenBank genomic record.

    Args:
        gb_path: Path to the GenBank file.

    Returns:
        Bio.SeqRecord.SeqRecord: The loaded GenBank record.

    Raises:
        FileNotFoundError: If the GenBank file is missing.
    """
    if not os.path.exists(gb_path):
        raise FileNotFoundError(f"Input file not found at: {gb_path}")

    record = SeqIO.read(gb_path, "genbank")
    return record


def extract_gene_nucleotide_sequences(
    record: SeqRecord,
) -> Dict[str, Dict[str, object]]:
    """Extract precise nucleotide sequence and genomic coordinates for every CDS.

    Utilizes SeqFeature.extract() to automatically handle forward/reverse strands
    and multipart spliced locations.

    Args:
        record: The parsed GenBank SeqRecord.

    Returns:
        A dictionary mapping gene symbols to their location coordinates,
        strand, and extracted DNA sequence.
    """
    gene_catalog: Dict[str, Dict[str, object]] = {}

    for feature in record.features:
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]

            # feature.extract(record.seq) automatically extracts spliced/stand-specific DNA!
            extracted_dna: Seq = feature.extract(record.seq)

            # Store coordinates as 1-based human-readable genomic positions
            start_pos = int(feature.location.start) + 1
            end_pos = int(feature.location.end)
            strand_symbol = "+" if feature.location.strand == 1 else "-"

            gene_catalog[gene_name] = {
                "start": start_pos,
                "end": end_pos,
                "strand": strand_symbol,
                "length_nt": len(extracted_dna),
                "dna_seq": str(extracted_dna),
                "is_compound": feature.location.parts is not None
                and len(feature.location.parts) > 1,
            }

    return gene_catalog


def check_viral_integration_break(
    gene_catalog: Dict[str, Dict[str, object]],
    breakpoint_coord: int,
) -> Dict[str, object]:
    """Analyze which viral gene is disrupted by a host-virus genomic integration breakpoint.

    In HPV-16, host genomic integration frequently disrupts the E2 regulatory gene
    (coordinates ~2756-3853 bp), which relieves repression on E6/E7 oncogenes.

    Args:
        gene_catalog: Extracted genes dictionary from extract_gene_nucleotide_sequences.
        breakpoint_coord: Genomic nucleotide coordinate (1-based) where the break occurs.

    Returns:
        A dictionary containing the disruption status, disrupted gene, and clinical impact.
    """
    disrupted_genes: List[str] = []

    for gene, info in gene_catalog.items():
        start = int(info["start"])
        end = int(info["end"])
        if start <= breakpoint_coord <= end:
            disrupted_genes.append(gene)

    if disrupted_genes:
        is_e2_disrupted = "E2" in disrupted_genes
        return {
            "breakpoint": breakpoint_coord,
            "disrupted": True,
            "genes": disrupted_genes,
            "clinical_consequence": (
                "Loss of E2 repression -> Hyperactivation of E6/E7 oncogenes (High Carcinogenic Risk)"
                if is_e2_disrupted
                else f"Interruption of {', '.join(disrupted_genes)}"
            ),
        }

    return {
        "breakpoint": breakpoint_coord,
        "disrupted": False,
        "genes": [],
        "clinical_consequence": "Breakpoint located in non-coding/intergenic region.",
    }


def display_integration_report(
    catalog: Dict[str, Dict[str, object]],
    break_result: Dict[str, object],
) -> None:
    """Print an academic summary of extracted gene coordinates and integration event.

    Args:
        catalog: Gene features catalog.
        break_result: Disruption analysis outcome.
    """
    print("=" * 70)
    print("      GENOMIC FEATURE COORDINATE & EXTRACTION AUDIT")
    print("=" * 70)
    print(
        f"{'Gene':<8} | {'Coordinates (bp)':<18} | {'Strand':<6} | {'Length (nt)':<12} | {'Spliced'}"
    )
    print("-" * 70)

    for gene, data in catalog.items():
        coord_str = f"{data['start']}..{data['end']}"
        spliced_str = "Yes (Compound)" if data["is_compound"] else "No"
        print(
            f"{gene:<8} | {coord_str:<18} | {str(data['strand']):<6} | {str(data['length_nt']):>8} nt  | {spliced_str}"
        )

    print("\n" + "=" * 70)
    print("        VIRAL HOST INTEGRATION DISRUPTION REPORT")
    print("=" * 70)
    print(f"Simulated Integration Breakpoint: {break_result['breakpoint']} bp")
    print(f"Gene Disrupted: {break_result['genes'] if break_result['disrupted'] else 'None'}")
    print(f"Biological Consequence:\n  -> {break_result['clinical_consequence']}")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    GB_FILE = "hpv16_genome.gb"
    BREAKPOINT = 3150  # Classical E2 hinge region integration site

    try:
        hpv_record = load_genomic_record(GB_FILE)
        genes_data = extract_gene_nucleotide_sequences(hpv_record)
        integration_status = check_viral_integration_break(
            genes_data, BREAKPOINT
        )
        display_integration_report(genes_data, integration_status)
    except Exception as err:
        print(f"[ERROR] Execution failed: {err}")
