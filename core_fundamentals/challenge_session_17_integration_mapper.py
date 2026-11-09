"""Session 17 Challenge: HPV-16 Genomic Feature Integration Mapper.

This script parses the HPV-16 GenBank file, builds a coordinate map of CDS features,
and scans clinical breakpoint locations to evaluate gene disruption risks (ZJU/PubMed).
"""

from typing import Dict, List, Tuple
from Bio import SeqIO


def build_cds_coordinate_map(gb_file: str) -> Dict[str, Tuple[int, int]]:
    """Parse GenBank record and map each gene name to its (start, end) coordinates."""
    record = SeqIO.read(gb_file, "genbank")
    gene_coordinate_map: Dict[str, Tuple[int, int]] = {}

    for feature in record.features:
        if feature.type == "CDS":
            gene_name = feature.qualifiers.get(
                "gene", feature.qualifiers.get("protein_id", ["unknown"])
            )[0]
            start = int(feature.location.start) + 1  # Convert to 1-based
            end = int(feature.location.end)          # Biopython end is exclusive
            gene_coordinate_map[gene_name] = (start, end)

    if not gene_coordinate_map:
        raise ValueError("Generated coordinate map is empty.")

    return gene_coordinate_map


def scan_clinical_breakpoints(
    cds_map: Dict[str, Tuple[int, int]], breakpoints: List[int]
) -> List[Dict[str, str]]:
    """Scan clinical breakpoints against gene coordinates to assess disruption risks."""
    integration_results: List[Dict[str, str]] = []

    for bp in breakpoints:
        disrupted_genes = []
        for gene_name, (start, end) in cds_map.items():
            if start <= bp <= end:
                disrupted_genes.append(gene_name)

        if disrupted_genes:
            is_e2_disrupted = "E2" in disrupted_genes
            status = (
                "HIGH_RISK (E2 Disrupted)"
                if is_e2_disrupted
                else f"MODERATE_RISK ({', '.join(disrupted_genes)} disrupted)"
            )
            integration_results.append(
                {"breakpoint": str(bp), "status": status}
            )
        else:
            integration_results.append(
                {"breakpoint": str(bp), "status": "LOW_RISK (Intergenic)"}
            )

    return integration_results


def report(results: List[Dict[str, str]]) -> None:
    """Print an academic formatted report of clinical breakpoint scans."""
    print("=" * 65)
    print("         HPV-16 CLINICAL BREAKPOINT INTEGRATION REPORT")
    print("=" * 65)
    print(f"{'Breakpoint (bp)':<18} | {'Clinical Risk & Status'}")
    print("-" * 65)

    for item in results:
        bp = item["breakpoint"]
        status = item["status"]
        print(f"{bp:<18} | {status}")

    print("=" * 65 + "\n")


if __name__ == "__main__":
    GB_FILE = "hpv16_genome.gb"
    BREAKPOINTS = [450, 860, 2700, 3150, 4200, 7100]

    try:
        cds_map = build_cds_coordinate_map(GB_FILE)
        scan_results = scan_clinical_breakpoints(cds_map, BREAKPOINTS)
        report(scan_results)
    except Exception as err:
        print(f"[ERROR] Pipeline execution failed: {err}")
