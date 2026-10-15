"""
Session 07: Physicochemical and Compositional Analysis of Proteins
Author: Surena Ghadam
Description: Analyzes amino acid frequencies and biochemical groups of translated proteins (e.g., Human Insulin).
"""
from session_05_fasta_parser import read_fasta 
from def_dna_to_protein import translate_dna

def analyze_amino_acids(protein_seq):
    aa_count={}
    
    pro_set=set(protein_seq)
    for aa in pro_set :
        aa_valiu=protein_seq.count(aa)
        aa_percent=( aa_valiu/len(protein_seq))*100
        aa_count[aa]=aa_valiu,aa_percent 
    return aa_count
def analyze_amino_acids_type(protein_seq):
    
    aa_type={}
    Hydrophobic_count=0
    Polar_Uncharged_count=0
    Acidic_Negative_Charge_count=0
    Basic_Positive_Charge_count=0

    Hydrophobic_aa={"A", "V", "I", "L", "M", "F", "W", "P"}
    Polar_Uncharged= {"S", "T", "C","N", "Q", "Y","G"}
    Acidic_Negative_Charge={"D", "E"}
    Basic_Positive_Charge={"R", "H", "K"}

    for aa in protein_seq :
        
        
        if aa in Hydrophobic_aa :
            Hydrophobic_count+=1
        elif aa in Polar_Uncharged :
            Polar_Uncharged_count+=1
        elif aa in Acidic_Negative_Charge :
            Acidic_Negative_Charge_count+=1
        elif aa in Basic_Positive_Charge:
            Basic_Positive_Charge_count+=1
    
    
    aa_Hydrophobic_percent=( Hydrophobic_count/len(protein_seq))*100
    aa_Polar_Uncharged_percent=( Polar_Uncharged_count/len(protein_seq))*100
    aa_Acidic_Negative_Charge_percent=( Acidic_Negative_Charge_count/len(protein_seq))*100
    aa_Basic_Positive_Charge_percent=( Basic_Positive_Charge_count/len(protein_seq))*100

    

    aa_type["Hydrophobic (A, V, I, L, M, F, W, P)"]=Hydrophobic_count,aa_Hydrophobic_percent
    aa_type["Polar Uncharged (S, T, C, N, Q, Y, G)"]=Polar_Uncharged_count,aa_Polar_Uncharged_percent
    aa_type["Acidic / Negative Charge (D, E)"]=Acidic_Negative_Charge_count,aa_Acidic_Negative_Charge_percent
    aa_type["Basic / Positive Charge (R, H, K)"]=Basic_Positive_Charge_count,aa_Basic_Positive_Charge_percent
    
    return aa_type




DNA_anti_fasta=read_fasta("human_insulin.fasta")
protein_insulin=translate_dna(DNA_anti_fasta)
analysis_result = analyze_amino_acids(protein_insulin)
analysis_type_result=analyze_amino_acids_type(protein_insulin)

print("\n[1] Amino Acid Composition (Count & Percentage):")
for aa , (Count,pct) in analysis_result.items() :
    print(f"- {aa} : {Count} ({pct:.2f})")
print("\n[2] Biochemical Classification:")
for group,(Count,pct) in analysis_type_result.items():
    print(f" {group} : {Count} , ({pct:.2f})")
