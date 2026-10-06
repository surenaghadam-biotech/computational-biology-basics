"""
Session 02: Hemoglobin Subunit Beta (HBB) Point Mutation Simulation
Topic: Python Lists, Mutability, and Sickle Cell Anemia (HbS) SNP Modeling
Author: Surena Ghadam

Description:
    Models the canonical missense point mutation in the HBB gene:
    A -> T transversion at codon 6 (GAG -> GTG), leading to the
    E6V (Glu -> Val) substitution characteristic of Sickle Cell Anemia.
"""

# Wild-type sequence: Codons 1 through 8 of Human HBB
# Codons: 1:ATG | 2:GTG | 3:CAC | 4:CTG | 5:ACT | 6:CCT | 7:GAG (Target) | 8:GAG
hbb_wildtype = "ATGGTGCAACTGACTCCTGAGGAG"

# 1. Convert immutable string to mutable list
list_hbb_wildtype = list(hbb_wildtype)

# 2. Introduce the missense mutation (A -> T) at codon 7 position 2 (Index -5 / Index 19)
list_hbb_wildtype[-5] = 'T'

# 3. Reconstruct the mutated DNA string
hbb_sickle = "".join(list_hbb_wildtype)

# 4. Extract target codon (HbA vs HbS) for biological validation
codon_wildtype = hbb_wildtype[-6:-3]
codon_sickle = hbb_sickle[-6:-3]

# 5. Formatted terminal output
print("=" * 65)
print("🧬 HBB MUTATION ANALYSIS: SICKLE CELL ANEMIA (HbS)")
print("=" * 65)
print(f"Wild-Type DNA (HbA) : {hbb_wildtype}")
print(f"Mutant DNA    (HbS) : {hbb_sickle}")
print("-" * 65)
print(f"Target Codon Before : {codon_wildtype} -> Glutamic Acid (Glu)")
print(f"Target Codon After  : {codon_sickle} -> Valine (Val) [Pathogenic SNP]")
print("=" * 65)
