"""
Session 10: Automated Genomic Data Retrieval via NCBI Entrez API
Author: Surena Ghadam
Project: Computational Biology Basics

Objectives:
1. Configure Entrez access parameters (Email compliance).
2. Fetch wild-type Human Hemoglobin Subunit Beta (HBB) mRNA (NM_000518.5).
3. Parse GenBank/FASTA directly from the network stream into SeqRecord.
4. Export the downloaded sequence locally as a standard FASTA file.
"""

from Bio import Entrez
from Bio import SeqIO

# ------------------------------------------------------------------------------
# 1. API Configuration & Authentication
# ------------------------------------------------------------------------------
# NCBI requires a contact email to track query rates and contact developers.
Entrez.email = "surenaghadam2007@gmail.com"  # Replace or keep standard email

# Target accession: NM_000518.5 (Homo sapiens hemoglobin subunit beta mRNA)
accession_id = "NM_000518.5"
database = "nucleotide"

print(f"Connecting to NCBI database [{database}]...")
print(f"Requesting Accession ID: {accession_id}")

# ------------------------------------------------------------------------------
# 2. Network Stream & Fetching (efetch)
# ------------------------------------------------------------------------------
# rettype="fasta": return in FASTA format (can also be "gb" for GenBank)
# retmode="text": plain text stream
with Entrez.efetch(db=database, id=accession_id, rettype="fasta", retmode="text") as handle:
    # SeqIO.read parses the network buffer directly without temporary disk files
    record = SeqIO.read(handle, "fasta")

# ------------------------------------------------------------------------------
# 3. Inspecting the Retrieved Sequence
# ------------------------------------------------------------------------------
print("\n--- Download Successful ---")
print(f"Record ID:          {record.id}")
print(f"Description:        {record.description}")
print(f"Sequence Length:    {len(record.seq)} bp")
print(f"First 30 bases:     {record.seq[:30]}...")

# ------------------------------------------------------------------------------
# 4. Caching Data to Local Disk
# ------------------------------------------------------------------------------
output_filename = "human_hbb_mrna.fasta"
SeqIO.write(record, output_filename, "fasta")
print(f"\nSequence cached locally as: {output_filename}")
