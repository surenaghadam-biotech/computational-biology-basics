# Computational Biology & Sequence Analytics Research Suite

[![Python](https://img.shields.io/badge/Python-3.8%2B-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![Biopython](https://img.shields.io/badge/Biopython-1.83-00599C.svg)](https://biopython.org/)
[![NCBI Entrez](https://img.shields.io/badge/NCBI-Entrez%20API-blue.svg)](https://www.ncbi.nlm.nih.gov/home/develop/api/)
[![Code Style](https://img.shields.io/badge/Code%20Style-PEP%208-orange.svg)](https://peps.python.org/pep-0008/)

An open-source computational genomics repository dedicated to reproducing published molecular biology literature, parsing NCBI databases, and developing in-silico pipelines for sequence alignment and physicochemical feature extraction.

---

## 🏛️ Grounded Literature & Case Studies

This repository implements and reproduces methodologies from peer-reviewed biomedical literature:

### 1. Comparative Insulin Genomics (*Sus scrofa* vs. *Homo sapiens*)
* **Reference Study:** *Evolution of the insulin molecule: insights into structure-activity relationships.*
* **Principal Investigator:** Dr. J. Michael Conlon (Department of Biochemistry, Faculty of Medicine and Health Sciences, UAE University).
* **Identifier:** [PMID: 11445250](https://pubmed.ncbi.nlm.nih.gov/11445250/) | *Peptides*, 2001.
* **In-Silico Reproduction:**
  * Implemented pairwise dot-plot matrix alignment between porcine and human preproinsulin.
  * Verified conservation across A-chain and B-chain domains and calculated sequence identity metrics.
  * Truncated 5' UTR sequences via codon-aware indexing to isolate exact biological CDS regions.

### 2. Hemoglobin Subunit Beta (*HBB*) Mutation Analysis & RefSeq Pipelines
* **Target RefSeq:** NCBI `NM_000518.5` (*Homo sapiens* hemoglobin subunit beta mRNA).
* **Pathogenic Variant Modeling:** Simulates the sickle cell single nucleotide polymorphism (classical `E6V` / HGVS `p.Glu7Val`, coding DNA: `c.20A>T`, mRNA: `GAG -> GUG`) and maps downstream physicochemical consequences on the $\beta$-globin chain.

* **Computational Focus:**
  * Automated streaming of complete nucleotide and GenBank flat files via `Bio.Entrez`.
  * Open Reading Frame (ORF) scanning and reverse complement validation ($5' \to 3'$).

---

## 📁 Repository Architecture
```text
├── core_fundamentals/                   # Progressive algorithmic sessions (PEP 8 compliant)
│   ├── session_01_dna_basics_en.py      # DNA transcription and base quantification
│   ├── session_02_hbb_sickle_cell_en.py # In-silico modeling of sickle cell mutation
│   ├── session_03_tp53_scanner.py       # TP53 tumor suppressor ORF scanner
│   ├── session_07_protein_analytics.py  # Physicochemical properties and AA frequency
│   ├── session_09_biopython_intro.py    # Biopython Bio.Seq and gc_fraction integration
│   └── session_10_ncbi_entrez.py        # Automated NCBI RefSeq retrieval via network streams
│
├── projects/
│   ├── 01_insulin_comparative_genomics/ # Case study reproducing PMID: 11445250
│   │   ├── data/                        # Raw & processed FASTA sequences
│   │   ├── scripts/                     # Dot-plot generation and identity matrices
│   │   └── results/                     # Visual sequence homology matrices
│   │
│   └── 02_hbb_mutation_and_entrez/      # Case study on Human HBB transcript
│       ├── data/                        # RefSeq records and CDS extractions
│       └── scripts/                     # Challenge workflows and translation validation
│
├── .gitignore
├── requirements.txt
└── README.md

```
## 🧬 Core Mathematical Formulations

  ### GC-Content Metric & Hybridization Stability Proxy
  The global nucleotide composition—used as an empirical proxy for duplex thermal stability and melting temperature ($T_m$)—is defined as:


$$\text{GC-Content (\%)} = \left( \frac{n_G + n_C}{L} \right) \times 100$$

where $L$ denotes the total sequence length:
$$L = n_A + n_T + n_U + n_G + n_C$$

### Pairwise Dot-Plot Homology Matrix
Given sequences $S_1$ of length $M$ and $S_2$ of length $N$, the binary pairwise alignment matrix $D \in \{0, 1\}^{M \times N}$ is computed over indices $1 \le i \le M$ and $1 \le j \le N$ as:

$$D_{i,j} = \begin{cases} 1 & \text{if } S_1[i] = S_2[j] \\ 0 & \text{otherwise} \end{cases}$$

---

## 🚀 Environment Setup & Reproducibility

```bash
# Clone the repository
git clone https://github.com/surenaghadam-biotech/computational-biology-basics.git
cd computational-biology-basics

# Install dependencies
pip install -r requirements.txt
```

---

## 👨💻 Author & Research Inquiries
**Surena Ghadam**  
*Aspiring Computational Biologist & Bioinformatics Researcher*  
* GitHub: [@surenaghadam-biotech](https://github.com/surenaghadam-biotech)


