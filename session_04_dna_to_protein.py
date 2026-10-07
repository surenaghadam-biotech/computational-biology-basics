"""
Session 04: DNA to Protein Translation Engine
Topic: Python Dictionaries, Key-Value Lookup, and Codon Translation
Target: Human HBB and TP53 Coding Sequences
Author: Surena Ghadam
"""

# The Standard Genetic Code Dictionary (DNA Codons -> Single-Letter Amino Acids)
# '*' denotes a STOP codon.
GENETIC_CODE = {
    # Phenylalanine & Leucine
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    # Leucine
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine & Start (Methionine)
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    # Proline
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine & Stop Codons
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    # Histidine & Glutamine
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    # Asparagine & Lysine
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    # Aspartate & Glutamate
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    # Cysteine, Stop Codon & Tryptophan
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    # Arginine
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    # Serine & Arginine
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    # Glycine
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}
tp53_wd = "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGA"
tp53_list=list(tp53_wd)
tp53_list[523]='A'
tp53_mut="".join(tp53_list)
count_codon=0
protein_wd = ""
protein_mut = ""
for i in range(0,len(tp53_wd),3):
    codon_wd=tp53_wd[i:i+3]
    codon_mut=tp53_mut[i:i+3]
    count_codon=count_codon+1
    codon_wd_aa=GENETIC_CODE[codon_wd]
    codon_mut_aa=GENETIC_CODE[codon_mut]
    if codon_wd_aa=="*":
        print("STOP TRANSLATE")
        break
    protein_wd += codon_wd_aa
    protein_mut += codon_mut_aa
    if codon_wd!=codon_mut:
        print(f"Warning! from codon{count_codon}, we have: {codon_wd_aa}->{codon_mut_aa} ")
print(f"pro_wt : {protein_wd} with: {len(protein_wd)} aa | pro_mut : {protein_mut} with: {len(protein_mut)} aa")

aa175_wt = protein_wd[174]
aa175_mut = protein_mut[174]

print(f"aa-175 : {aa175_wt} at WT and {aa175_mut} at Mutant")
