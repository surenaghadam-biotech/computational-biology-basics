"""
Session 02: Hemoglobin Subunit Beta (HBB) Point Mutation Simulation
Topic: Python Lists, Mutability, and Sickle Cell Anemia (HbS) SNP Modeling
Author: Surena Ghadam

Description:
    Models the canonical missense point mutation in the HBB gene:
    A -> T transversion at codon 6 (GAG -> GTG), leading to the
    E6V (Glu -> Val) substitution characteristic of Sickle Cell Anemia.
"""
hbb_wildtype = "ATGGTGCAACTGACTCCTGAGGAG"

list_hbb_wildtype=list(hbb_wildtype)
list_hbb_wildtype[-5]='T'
hbb_sickle="".join(list_hbb_wildtype)
print(f"توالی قبل از جهش :{hbb_wildtype} | توالی پس از جهش : {hbb_sickle}")