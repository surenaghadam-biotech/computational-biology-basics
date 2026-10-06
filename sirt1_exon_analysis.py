"""
SIRT1 Exon Characterization & Translation Initiation Check

Author: Surena Ghadam
Context: Longevity Genomics / Sirtuin Pathway Analysis

Description:
    Parses a target exon of the SIRT1 gene to validate the canonical
    ATG translation initiation codon and compute nucleotide compositional
    metrics (AT/GC ratio).
"""
sirt1_exon = "ATGGCGGACGAGGCGGCGCTCGCCCTGCAGCCCGGCGGTTCCAGCTCGGCTCGG"

sirt1_len=len(sirt1_exon)
print(f'طول رشته کد کننده:{sirt1_len}نوکلیوتید')

start_gene=sirt1_exon[0:3]
end_gene=sirt1_exon[-3: ]
print(f'کد اغاز ژن:{start_gene} و کد پایان ژن: {end_gene}  ')

if start_gene=='ATG' :
    print('کد اغاز برای امینو اسید متیونین است')
else:
    print('این  کد اغاز برای ترجمه نیست نه')

count_a=sirt1_exon.count('A')
count_t=sirt1_exon.count('T')
count_c=sirt1_exon.count('C')
count_g=sirt1_exon.count('G')
AT_gene= count_a + count_t
GC_gene= count_g + count_c

at_gc_ratio = AT_gene / GC_gene
print(f"نسبت AT به GC: {at_gc_ratio:.2f}")


at_percent = (AT_gene / sirt1_len) * 100
gc_percent = (GC_gene / sirt1_len) * 100
print(f"محتوای AT: {at_percent:.1f}% | محتوای GC: {gc_percent:.1f}%")


