"""
Session 03: Genome Scanner and Reading Frame Analysis
Topic: For Loops, Range Stepping, and Stop Codon Identification
Target: TP53 Tumor Suppressor Exon Segment
Author: Surena Ghadam
"""
tp53_wd = "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCTGTCCCCCTTGCCGTCCCAAGCAATGGATGATTTGATGCTGTCCCCGGACGATATTGAACAATGGTTCACTGAAGACCCAGGTCCAGATGAAGCTCCCAGAATGCCAGAGGCTGCTCCCCCCGTGGCCCCTGCACCAGCAGCTCCTACACCGGCGGCCCCTGCACCAGCCCCCTCCTGGCCCCTGTCATCTTCTGTCCCTTCCCAGAAAACCTACCAGGGCAGCTACGGTTTCCGTCTGGGCTTCTTGCATTCTGGGACAGCCAAGTCTGTGACTTGCACGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGGTCTGGCCCCTCCTCAGCATCTTATCCGAGTGGAAGGAAATTTGCGTGTGGAGTATTTGGATGACAGAAACACTTTTCGACATAGTGTGGTGGTGCCCTATGAGCCGCCTGAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAAGACTCCAGTGGTAATCTACTGGGACGGAACAGCTTTGAGGTGCGTGTTTGTGCCTGTCCTGGGAGAGACCGGCGCACAGAGGAAGAGAATCTCCGCAAGAAAGGGGAGCCTCACCACGAGCTGCCCCCAGGGAGCACTAAGCGAGCACTGCCCAACAACACCAGCTCCTCTCCCCAGCCAAAGAAGAAACCACTGGATGGAGAATATTTCACCCTTCAGATCCGTGGGCGTGAGCGCTTCGAGATGTTCCGAGAGCTGAATGAGGCCTTGGAACTCAAGGATGCCCAGGCTGGGAAGGAGCCAGGGGGGAGCAGGGCTCACTCCAGCCACCTGAAGTCCAAAAAGGGTCAGTCTACCTCCCGCCATAAAAAACTCATGTTCAAGACAGAAGGGCCTGACTCAGACTGA"
tp53_list=list(tp53_wd)
tp53_list[523]='A'
tp53_mut="".join(tp53_list)

total_codon=0
mutation_count = 0

for i in range(0,len(tp53_wd),3):
    codon_wd=tp53_wd[i:i+3]
    codon_mut=tp53_mut[i:i+3]
    codon_num=(i//3)+1 
    total_codon=total_codon+1
    if codon_wd != codon_mut :
        print(f'🚨 MUTATION at Codon {codon_num} (nt {i+1}-{i+3}): {codon_wd} -> {codon_mut} [R175H HOTSPOT]')
        mutation_count = mutation_count + 1
    if codon_mut=="ATG" and total_codon==1  :
        print("[-----START CODON-----]")
    if codon_mut in ["TAA", "TAG", "TGA"]:
         if total_codon == 394:
             print(f"Codon {codon_num} | Pos: {i+1}-{i+3} | Sequence: {codon_mut} | Status: [NORMAL STOP CODON]")
         else:
             print(f"🚨 PREMATURE STOP CODON at Codon {codon_num} (nt {i+1}-{i+3}): {codon_mut}")


# Report Section (English)
print("-" * 40)
print(f"Gene Length : {len(tp53_wd)} bp | Total Codons: {total_codon}")

# Wild-Type GC Calculation
count_c = tp53_wd.count("C")
count_g = tp53_wd.count("G")
gc_sum = count_c + count_g
percent_gc = (gc_sum / len(tp53_wd)) * 100
print(f"Wild-Type GC Content : {percent_gc:.2f}%")

# Mutant GC Calculation
count_c2 = tp53_mut.count("C")
count_g2 = tp53_mut.count("G")
gc_sum2 = count_c2 + count_g2
percent_gc2 = (gc_sum2 / len(tp53_mut)) * 100
print(f"Mutant GC Content    : {percent_gc2:.2f}%")

print(f"Total Mutations      : {mutation_count}")
print("-" * 40)




        
