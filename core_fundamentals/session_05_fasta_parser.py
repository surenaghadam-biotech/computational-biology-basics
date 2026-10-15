"""
Session 05: FASTA File Parser and mRNA Sequence Extractor
Topic: File I/O, Line Stripping, Header Filtering, and Reusable Functions
Target: Human Insulin (INS) mRNA (NCBI: NM_000207.3)
Author: Surena Ghadam
"""


def read_fasta(file_path):
    """Reads a FASTA file and extracts the pure nucleotide sequence.

    Skips metadata headers starting with '>'.
    """
    dna = ""
    with open(file_path, "r") as file:
        for line in file:
            clean_line = line.strip()
            if not clean_line.startswith(">"):
                dna += clean_line
    return dna


# --- Main Execution ---
file_name = "human_insulin.fasta"
insulin_seq = read_fasta(file_name)

print("=" * 45)
print("     HUMAN INSULIN (INS) FASTA REPORT")
print("=" * 45)
print(f"Source File          : {file_name}")
print(f"Total Sequence Length: {len(insulin_seq)} bp")
print(f"First 50 Nucleotides : {insulin_seq[:50]}")
print(f"Last 20 Nucleotides  : {insulin_seq[-20:]}")
print("=" * 45)



