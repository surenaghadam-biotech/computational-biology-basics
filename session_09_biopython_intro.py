from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

# 1. Read FASTA file and wrap it into a SeqRecord object using SeqIO.read
record = SeqIO.read("pig_insulin.fasta", "fasta")

# Display record metadata and sequence properties
print("Gene ID     :", record.id)
print("Header Desc :", record.description)
print("Seq Length  :", len(record.seq))
print("First 10 bp :", record.seq[:10])
print("Seq Type    :", type(record.seq))

# Store sequence object in a variable for convenience
my_dna = record.seq

print("\n--- Biopython Biological Methods Test ---")

# 1. Original sequence slice (first 10 bp)
print("Original DNA (first 10 bp):", my_dna[:10])

# 2. Complement strand
print("Complement strand         :", my_dna[:10].complement())

# 3. Reverse Complement strand
print("Reverse Complement        :", my_dna[:10].reverse_complement())

# 4. Transcription to RNA
my_rna = my_dna[:10].transcribe()
print("Transcribed to RNA        :", my_rna)

print("\n--- Translation Test ---")

# Translate full pig mRNA sequence to amino acids
my_protein = my_dna.translate()
print("Full pig protein:")
print(my_protein)

# Translate up to the first stop codon (to_stop=True)
my_protein_to_stop = my_dna.translate(to_stop=True)
print("\nProtein up to stop codon:")
print(my_protein_to_stop)
print("Resulting protein length:", len(my_protein_to_stop), "amino acids")

print("\n--- Sequence Analytics with Bio.SeqUtils ---")

# Calculate GC content percentage for the pig insulin gene
gc_value = gc_fraction(my_dna) * 100
print(f"Pig insulin GC content: {gc_value:.2f}%")

# Count exact nucleotide frequencies using count() method
count_a = my_dna.count("A")
count_t = my_dna.count("T")
count_c = my_dna.count("C")
count_g = my_dna.count("G")

print("\nExact nucleotide counts:")
print(f"A: {count_a} | T: {count_t} | C: {count_c} | G: {count_g}")

print("\n--- Precise Biological Translation (From start codon ATG) ---")

# 1. Find the index of the first start codon (ATG)
start_index = my_dna.find("ATG")
print(f"First start codon (ATG) found at index: {start_index}")

# 2. Trim sequence from the start codon onward (removing 5' UTR)
coding_dna = my_dna[start_index:]
print(f"Coding sequence length: {len(coding_dna)} nucleotides")

# 3. Translate trimmed sequence up to stop codon
real_insulin_protein = coding_dna.translate(to_stop=True)

print("\nReal pig insulin protein (starting with M):")
print(real_insulin_protein)
print(f"Exact protein length: {len(real_insulin_protein)} amino acids")

print("\n--- Export Final Protein to FASTA File ---")

# 1. Wrap protein sequence into a new SeqRecord object
protein_record = SeqRecord(
    real_insulin_protein,
    id="pig_insulin_protein",
    description="Sus scrofa translated preproinsulin (120 aa)",
)

# 2. Save to disk using SeqIO.write
output_filename = "pig_insulin_protein.fasta"
SeqIO.write(protein_record, output_filename, "fasta")

print(f"File successfully saved: {output_filename}")
