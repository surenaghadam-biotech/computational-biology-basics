"""
Analyze HBB mRNA Sequence and Extract Coding Sequence (CDS)
Author: Surena Ghadam
Project: Computational Biology Basics (02_hbb_mutation_and_entrez)

Pipeline Workflow:
1. Load human beta-globin (HBB, NM_000518.5) full mRNA sequence from FASTA.
2. Quantify basic genomic properties (sequence length and GC content).
3. Derive reverse complement cDNA for downstream PCR primer design.
4. Scan canonical start codon (ATG/AUG) and trim upstream 5' UTR.
5. Translate the open reading frame (ORF) to protein and export results.
"""

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


def run_pipeline() -> None:
    # -------------------------------------------------------------------------
    # Setup Dynamic Relative Paths
    # -------------------------------------------------------------------------
    base_dir = Path(__file__).resolve().parent.parent
    input_mrna_path = base_dir / "data" / "human_hbb_mrna.fasta"
    results_dir = base_dir / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    output_protein_path = results_dir / "hbb_protein.fasta"

    print("=" * 60)
    print("HBB mRNA ANALYSIS & CDS EXTRACTION PIPELINE")
    print("=" * 60)

    # -------------------------------------------------------------------------
    # Step 1: Ingest Sequence Record from FASTA
    # -------------------------------------------------------------------------
    if not input_mrna_path.exists():
        raise FileNotFoundError(f"Input file not found at: {input_mrna_path}")

    record = SeqIO.read(input_mrna_path, "fasta")
    mrna_seq = record.seq
    print(f"[1] Loaded Record ID: {record.id}")
    print(f"    Description: {record.description}")

    # -------------------------------------------------------------------------
    # Step 2: Genomic Analytics (Length & GC Content)
    # -------------------------------------------------------------------------
    seq_length = len(mrna_seq)
    gc_percent = gc_fraction(mrna_seq) * 100
    print(f"[2] Sequence Length: {seq_length} nt")
    print(f"    GC Content: {gc_percent:.2f}%")

    # -------------------------------------------------------------------------
    # Step 3: Reverse Complement cDNA Generation
    # -------------------------------------------------------------------------
    rev_comp = mrna_seq.reverse_complement()
    print(f"[3] Reverse Complement (3' end preview): ...{rev_comp[-30:]}")

    # -------------------------------------------------------------------------
    # Step 4: Scan Start Codon & Trim 5' UTR
    # -------------------------------------------------------------------------
    start_idx = mrna_seq.find("ATG")
    if start_idx == -1:
        start_idx = mrna_seq.find("AUG")

    if start_idx == -1:
        raise ValueError("Error: Canonical start codon (ATG/AUG) not found.")

    cds_candidate = mrna_seq[start_idx:]
    print(f"[4] Canonical start codon identified at position: {start_idx} (0-based index)")
    print(f"    5' UTR length trimmed: {start_idx} nt")

    # Trim any trailing partial codon nucleotides to ensure exact triplets
    valid_len = len(cds_candidate) - (len(cds_candidate) % 3)
    clean_cds = cds_candidate[:valid_len]

    # -------------------------------------------------------------------------
    # Step 5: Translate to Protein & Export Record
    # -------------------------------------------------------------------------
    protein_seq = clean_cds.translate(to_stop=True)
    print(f"[5] Translated Amino Acid Chain Length: {len(protein_seq)} aa")
    print(f"    Protein Sequence: {protein_seq}")

    protein_record = SeqRecord(
        seq=protein_seq,
        id="HBB_HUMAN_TRANSLATED",
        description="Hemoglobin subunit beta [Homo sapiens] translated CDS"
    )

    SeqIO.write(protein_record, output_protein_path, "fasta")
    print("-" * 60)
    print(f"[+] Output protein saved successfully to:\n    {output_protein_path}")
    print("=" * 60)


if __name__ == "__main__":
    run_pipeline()
