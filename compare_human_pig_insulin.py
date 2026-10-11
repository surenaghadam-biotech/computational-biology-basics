"""
Research Project: Comparative Insulin Sequence Analysis (Human vs. Pig)
PubMed Reference: PMID: 11445250 (Evolution of the insulin molecule)
Author: Surena Ghadam
Description: Compares human and pig insulin protein sequences, aligns them using 
             a dot plot matrix, and identifies evolutionary mutations (e.g., B30 position).
"""
from session_05_fasta_parser import read_fasta
from def_dna_to_protein import translate_dna
from def_dotplot_alignment import generate_dot_plot
from calculate_identity_score import calculate_identity_score

fasta_file_P="pig_insulin.fasta"
pig_mRNA=read_fasta(fasta_file_P)
pig_pro=translate_dna(pig_mRNA)
fasta_file_H ="human_insulin.fasta"
human_mRNA =read_fasta(fasta_file_H)
human_pro=translate_dna(human_mRNA)

generate_dot_plot(pig_pro,human_pro)
score=calculate_identity_score(pig_pro,human_pro)

print(f"Identity Score: {score:.2f}%")




