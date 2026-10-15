"""
Challenge Session 09: HBB mRNA Sequence Analysis and CDS Extraction
Author: Surena Ghadam
Project: Computational Biology Basics

Pipeline Workflow:
1. Load raw human beta-globin (HBB) partial mRNA sequence.
2. Quantify basic genomic properties (length and GC content).
3. Derive reverse complement cDNA for potential PCR primer design.
4. Scan and trim 5' UTR up to the canonical start codon (AUG).
5. Translate the open reading frame (ORF) to protein and export as FASTA.
"""

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# ==============================================================================
# Task 1: Genomic Analysis (Length & GC Content)
# ==============================================================================
raw_mrna_str = (
    "ACAUUUGCUUCUGACACAACUGUGUUCACUAGCAACCUCAAACAGACACC"
    "AUGGUGCACCUGACUCCUGAGGAGAAGUCUGCCGUUACUGCCCUGUGGGCAAGUGACUGAUCA"
)
raw_seq = Seq(raw_mrna_str)

# Calculate nucleotide length and GC percentage
raw_seq_len = len(raw_seq)
gc_percent = gc_fraction(raw_seq) * 100

print(f"Raw sequence length: {raw_seq_len} nt")
print(f"GC content: {gc_percent:.2f}%")

# ==============================================================================
# Task 2: Reverse Complement Generation (cDNA Strand)
# ==============================================================================
# Back-transcribe RNA to DNA, then derive the reverse complement strand
dna_seq = raw_seq.back_transcribe()
rev_comp_dna = dna_seq.reverse_complement()

print(f"\nReverse complement (5'->3' DNA): {rev_comp_dna}")

# ==============================================================================
# Task 3 & 4: Start Codon Scanning, 5' UTR Trimming, and Translation
# ==============================================================================
# Locate canonical start codon (AUG) to establish the correct reading frame
start_idx = raw_seq.find("AUG")

if start_idx != -1:
    # Trim upstream 5' Untranslated Region (5' UTR)
    coding_rna = raw_seq[start_idx:]
    
    # Translate coding sequence directly to the first stop codon
    protein_seq = coding_rna.translate(to_stop=True)
    
    if len(protein_seq) > 0:
        print(f"\nStart codon index: {start_idx}")
        print(f"First amino acid: {protein_seq[0]}")
        print(f"Translated protein: {protein_seq} (Length: {len(protein_seq)} aa)")
        
        # ======================================================================
        # Task 5: Exporting Translated Protein to FASTA Format
        # ======================================================================
        protein_record = SeqRecord(
            seq=protein_seq,
            id="HBB_HUMAN_PARTIAL",
            description="Hemoglobin subunit beta translated CDS"
        )
        
        output_filename = "hbb_protein.fasta"
        SeqIO.write(protein_record, output_filename, "fasta")
        print(f"\nSuccessfully exported protein to: {output_filename}")
    else:
        print("\nWarning: Translation yielded an empty sequence.")
else:
    print("\nError: Start codon (AUG) was not identified in the sequence.")
