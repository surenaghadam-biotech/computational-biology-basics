"""Session 16: Advanced GenBank Parsing and Genomic Feature Extraction.

Target: Automated extraction of CDS, protein products, and coordinates
for HPV-16 genome (NC_001526.4).
Designed for CSC China 2027 Resume Portfolio (ZJU benchmark).
"""

import os
from typing import Dict, List
from Bio import Entrez, SeqIO
from Bio.SeqFeature import SeqFeature
from utils import get_ncbi_email
 

Entrez.email = get_ncbi_email()

def fetch_genbank_record(accession_id: str, output_file: str) -> None:
    """Fetches a full GenBank record from NCBI Nucleotide and saves it locally.

    Skips download if file already exists.
    """
    if os.path.exists(output_file):
        print(
            f"[INFO] '{output_file}' already exists. Skipping online fetch."
        )
        return

    print(
        f"[INFO] Fetching '{accession_id}' in GenBank format from NCBI..."
    )
    with Entrez.efetch(
        db="nuccore", id=accession_id, rettype="gb", retmode="text"
    ) as handle:
        data = handle.read()
        
    with open(output_file, "w", encoding="utf-8") as out_handle:
        out_handle.write(data)
    print(f"[SUCCESS] Saved GenBank file to '{output_file}'.")


def parse_and_summarize_genbank(file_path: str) -> List[Dict[str, object]]:
    """Parses a GenBank file, summarizes top-level metadata, and extracts

    detailed information for each Coding DNA Sequence (CDS).
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    record = SeqIO.read(file_path, "genbank")

    print("\n" + "=" * 60)
    print(f"GENOME METADATA SUMMARY:")
    print(f"  Accession ID : {record.id}")
    print(f"  Organism     : {record.annotations.get('organism', 'N/A')}")#type of organism
    print(f"  Total Length : {len(record.seq):,} bp")# :, -> enrter "," at any 3rd seq: 7321->7,321
    print(f"  Topology     : {record.annotations.get('topology', 'linear')}")#circle/line
    print(f"  Features     : {len(record.features)} items detected")
    print("=" * 60 + "\n")

    cds_catalog: List[Dict[str, object]] = []

    for feature in record.features:
        feature: SeqFeature
        if feature.type == "CDS":

            gene_name = feature.qualifiers.get("gene", ["Unknown"])[0]
            product_name = feature.qualifiers.get(
                "product", ["Hypothetical"]
            )[0]
            protein_id = feature.qualifiers.get("protein_id", ["N/A"])[0]

            start_pos = int(feature.location.start)
            end_pos = int(feature.location.end)
            strand = "+" if feature.location.strand == 1 else "-"
            length_nt = end_pos - start_pos

            translation = feature.qualifiers.get("translation", [""])[0]
            protein_len = len(translation) if translation else length_nt // 3

            cds_info = {
                "gene": gene_name,
                "product": product_name,
                "protein_id": protein_id,
                "start": start_pos,
                "end": end_pos,
                "strand": strand,
                "length_nt": length_nt,
                "length_aa": protein_len,
            }
            cds_catalog.append(cds_info)

    return cds_catalog


def display_cds_table(cds_catalog: List[Dict[str, object]]) -> None:
    """Prints an aligned, publication-grade table of CDS features."""
    header = (
        f"{'Gene':<8} | {'Strand':<6} | {'Start':<6} | {'End':<6} | "
        f"{'NT Len':<6} | {'AA Len':<6} | {'Protein ID':<12} | {'Product'}"
    )
    print(header)
    print("-" * len(header))

    for item in cds_catalog:
        row = (
            f"{str(item['gene']):<8} | "
            f"{str(item['strand']):<6} | "
            f"{str(item['start']):<6} | "
            f"{str(item['end']):<6} | "
            f"{str(item['length_nt']):<6} | "
            f"{str(item['length_aa']):<6} | "
            f"{str(item['protein_id']):<12} | "
            f"{str(item['product'])}"
        )
        print(row)
    print("-" * len(header))
    print(f"Total CDS Extracted: {len(cds_catalog)}\n")


if __name__ == "__main__":
    ACCESSION = "NC_001526.4"
    OUTPUT_FILE = "hpv16_genome.gb"

    fetch_genbank_record(ACCESSION, OUTPUT_FILE)

    extracted_cds = parse_and_summarize_genbank(OUTPUT_FILE)

    display_cds_table(extracted_cds)
