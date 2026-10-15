from session_05_fasta_parser import read_fasta
from def_dna_to_protein import translate_dna
fasta_file ="human_insulin.fasta"
insulin_mrna =read_fasta(fasta_file)
insulin_pro=translate_dna(insulin_mrna)
print(f"pro of insulin :{insulin_pro}")