# 🧬 In Silico Simulation of Aberrant Splicing in $\beta$-Thalassemia ($HBB$ IVS-II-654)

[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![Biopython 1.83](https://img.shields.io/badge/Biopython-1.83-brightgreen.svg)](https://biopython.org/)
[![NCBI GenBank](https://img.shields.io/badge/NCBI-U01317.1-orange.svg)](https://www.ncbi.nlm.nih.gov/nuccore/U01317.1)
[![Reference Paper](https://img.shields.io/badge/PubMed-PMC10425742-red.svg)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10425742/)

---

## 📌 Biological Background & Overview

$\beta$-Thalassemia is an inherited hemoglobinopathy characterized by absent or reduced synthesis of $\beta$-globin chains. One of the most prevalent pathogenic non-coding mutations in East Asia (specifically documented by researchers at **Zhejiang University, China**) is the deep-intronic point mutation:

$$\text{IVS-II-654}\ (C \rightarrow T)$$

### 🔬 Molecular Mechanism:
1. **Aberrant Splicing Activation:** The substitution introduces a de novo $5^\prime$ splice donor site at position 654 and simultaneously activates a cryptic $3^\prime$ splice acceptor site at nucleotide 579 of Intron 2.
2. **Pseudo-Exon Inclusion:** This leads to the insertion of a **73 bp aberrant cryptic exon** between Exon 2 and Exon 3 in the mature mRNA transcript.
3. **Frameshift & Premature Termination Codon (PTC):** The 73 bp insertion disrupts the natural open reading frame (ORF), generating a premature stop codon that truncates the normal 147-amino-acid chain to a non-functional 121-amino-acid peptide.
```text
Wild-Type Pre-mRNA:   [ Exon 1 ]------( Intron 1 )------[ Exon 2 ]-------------------( Intron 2 )-------------------[ Exon 3 ]

Mutant Pre-mRNA:      [ Exon 1 ]------( Intron 1 )------[ Exon 2 ]----[ Cryptic Exon (73 bp) ]----[ Exon 3 ]
(IVS-II-654 C>T)                                                           (579 to 652)
```
## 🎯 Project Objectives

- Dynamically stream the complete human $\beta$-globin gene cluster from NCBI GenBank (`U01317.1`) using Biopython's `Entrez` and `SeqIO`.
- Programmatically isolate and extract the 3 functional coding exons of the $HBB$ gene.
- Simulate the intronic single nucleotide variant (SNV) $654\ C \rightarrow T$ inside Intron 2.
- Characterize the exact 73 bp cryptic exon integration.
- Translate both Wild-Type and Mutant transcripts to demonstrate the molecular truncation via PTC.

---

## 📂 Project Structure
```text
03_zju_hbb_aberrant_splicing/
│
├── README.md                              # Scientific documentation and overview
└── scripts/
└── simulate_aberrant_splicing.py      # Core computational biology pipeline
```
## ⚙️ Computational Pipeline

The pipeline (`simulate_aberrant_splicing.py`) operates in five modular steps:

1. **NCBI Querying:** Queries GenBank using `Entrez.efetch` to retrieve metadata and nucleotide sequences.
2. **Feature Filtering:** Filters genomic features specifically matching `type == 'exon'` and `gene == ['HBB']`.
3. **In Silico Mutagenesis:** Introduces the $C \rightarrow T$ transition at position 654 in Intron 2.
4. **Cryptic Splicing Simulation:** Slices the 73 bp intronic fragment (`579:652`) and splices it between Exon 2 and Exon 3.
5. **Translation & Alignment:** Scans for the initiating `ATG` start codon and translates up to the first encountered stop codon (`to_stop=True`).

---

## 📊 Verification & Console Output

When executed, the script confirms biological concordance with empirical studies:
```text
we have 3 exons at this wd gene.
we have 2 intron at this wd gene
Length of cryptic exon: 73 bp

======================================================================
BIOINFORMATICS ANALYSIS REPORT: HBB IVS-II-654 SPLICING
======================================================================
 [+] Wild-Type Protein Length : 147 amino acids
 [+] Mutant Protein Length    : 121 amino acids (Truncated via PTC)
----------------------------------------------------------------------

[Wild-Type Protein Sequence]
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH

[Mutant Protein Sequence (Pathogenic)]
MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRGNNDTMYHASLHHSKE
======================================================================
```

### 🧬 Sequence Comparison Analysis:
- **Residues 1–104:** Identical between Wild-Type and Mutant (`MVHL...PENFR`).
- **Residues 105–121:** Distinct aberrant sequence (`GNNDTMYHASLHHSKE`) induced by the frameshift, terminating immediately at residue 121 due to a PTC.

---

## 🚀 How to Run

Ensure Python 3.8+ and Biopython are installed in your environment:
```bash
pip install biopython
```
Navigate to the repository root and execute the script:
```bash
python projects/03_zju_hbb_aberrant_splicing/scripts/simulate_aberrant_splicing.py
```
## 📚 References

- **Paper:** *Therapeutic genome editing of an aberrant splice site in $\beta$-thalassemia by CRISPR/Cas9 with multiple sgRNAs* (Yang F. et al., 2023). **Zhejiang University**. [PMC10425742](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10425742/).
- **NCBI Reference:** Human $\beta$-globin region on chromosome 11, GenBank Accession: [U01317.1](https://www.ncbi.nlm.nih.gov/nuccore/U01317.1).
