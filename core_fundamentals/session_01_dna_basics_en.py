
"""
Computational Bioinformatics Project — Session 1
Topic: Basic DNA sequence analysis, slicing, and GC content calculation
Author: Surena Ghadam
"""

# 1. Define variables for the target gene sequence
gene_name = "LONGEVITY_CANDIDATE_GENE_1"
dna_seq = "ATGCGATCGATCGCTAGCTAGCTAGCGCTAGCTAGCTAATCGATCGATCGATAG"

print("=" * 65)
print(f"🧬 Bioinformatics analysis for gene: {gene_name}")
print("=" * 65)

# 2. Get the sequence length 
seq_length = len(dna_seq)
print(f"1. Gene sequence length: {seq_length} base pairs (bp)")

# 3. Access a nucleotide directly 
first_base = dna_seq[0]  # The first nucleotide 
print(f"2. First base at the 5′ end: {first_base}")

# 4. Extract a subsequence 
start_codon = dna_seq[0:3]  
print(f"3. Extracted start codon: {start_codon}")

# 5. Read from the end of the sequence
stop_codon_candidate = dna_seq[-3:] 
print(f"4. Candidate stop codon at the 3′ end: {stop_codon_candidate}")

# 6. Count the occurrences of each nucleotide 
count_a = dna_seq.count("A")
count_t = dna_seq.count("T")
count_c = dna_seq.count("C")
count_g = dna_seq.count("G")

print("-" * 65)
print("5. Nucleotide counts:")
print(f"   - Adenine (A): {count_a}")
print(f"   - Thymine (T): {count_t}")
print(f"   - Cytosine (C): {count_c}")
print(f"   - Guanine (G): {count_g}")

# 7. Calculate GC content
gc_count = count_g + count_c
gc_content = (gc_count / seq_length) * 100

print("-" * 65)
print(f"6. Sequence GC content: {gc_content:.2f}%")

if gc_content > 50.0:
    print("   [Biological interpretation: GC-rich sequence; potentially higher structural stability and melting temperature]")
else:
    print("   [Biological interpretation: AT-rich sequence; may denature at a lower temperature]")

print("=" * 65)
