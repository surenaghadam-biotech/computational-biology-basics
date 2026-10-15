# 🧬 Comparative Genomics: Human vs. Sus scrofa (Porcine) Insulin

[![Topic: Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics%20%26%20Genomics-blue.svg)](#)
[![Python: 3.8+](https://img.shields.io/badge/Python-3.8%2B-brightgreen.svg)](#)
[![Status: Reproducible](https://img.shields.io/badge/Reproducibility-Verified-orange.svg)](#)

## 📌 Biological Overview
Prior to the advent of recombinant DNA technology and biosynthetic human insulin (*humulin*), porcine (*Sus scrofa*) insulin served as the therapeutic gold standard in clinical endocrinology. Despite evolutionary divergence, porcine insulin exhibits near-identical physiological activity in human metabolic pathways due to extreme structural conservation across the functional A and B chains.

This module implements an automated comparative genomics pipeline to quantify sequence homology, calculate exact Hamming distances, and generate a 2D pairwise identity dot-plot matrix.

---

## 🔬 Mathematical Formulations

### 1. Pairwise Identity Dot-Plot Matrix
Given sequences $S_1$ and $S_2$ of lengths $M$ and $N$ respectively, the binary identity matrix $D \in \{0, 1\}^{M \times N}$ is formulated as:

$$D_{i, j} = \begin{cases} 1 & \text{if } S_1[i] = S_2[j] \\ 0 & \text{if } S_1[i] \neq S_2[j] \end{cases}$$

Continuous diagonals along the identity axis ($i = j$) represent identical peptide segments, whereas horizontal or vertical shifts delineate insertion-deletion (indel) events.

### 2. Sequence Identity & Hamming Metric
For ungapped homologous sequences of equal length $L = 110\text{ aa}$, the global sequence identity is inversely related to the Hamming distance $d_H(S_1, S_2)$:

$$\text{Identity (\%)} = \left( 1 - \frac{d_H(S_1, S_2)}{L} \right) \times 100 = \left( \frac{1}{L} \sum_{k=1}^{L} D_{k, k} \right) \times 100$$

---

## 📊 Visual Analytics & Structural Interpretation

The generated pairwise dot-plot illustrates the alignment between the preproinsulin chains:

![Human vs Porcine Insulin Dot-Plot](dotplot_human_pig_insulin.png)

### Key Structural Insights:
- **High-Density Primary Diagonal:** A nearly continuous diagonal spans coordinates $(1, 1)$ to $(110, 110)$, interrupted by only 2 point variations, confirming extreme sequence conservation ($98.18\%$ identity).
- **Conserved Structural Scaffold:** Invariant cysteine residues govern intra- and inter-chain disulfide bridges ($A6-A11$, $A7-B7$, and $A20-B19$), safeguarding the core tertiary folding necessary for receptor recognition.
- **Sequence Divergence Points:** The two mismatches do not disrupt receptor affinity, notably preserving the bioactive A and B interfaces while exhibiting minor divergence (such as the C-terminal B30 residue: Threonine in humans vs. Alanine in pigs).


---

## 💻 Script Execution & Output

Execute the analysis pipeline within the project directory:
```bash
python calculate_identity_score.py
```
### Pipeline Terminal Output:
```text
[INFO] Loading sequences...
[INFO] Target: Homo sapiens (Human) Preproinsulin (Length: 110 aa)
[INFO] Target: Sus scrofa (Porcine) Preproinsulin (Length: 110 aa)
[COMPUTATION] Hamming Distance: 2 mismatches
[COMPUTATION] Global Sequence Identity: 98.18%
[OUTPUT] Binary dot-plot rendered: dotplot_human_pig_insulin.png
```

## 📚 References & Data Sources
NCBI Entrez Database: Human Preproinsulin (INS, Accession: NP_000198.1) & Porcine Insulin (INS, Accession: NP_999424.1).
Conlon, J. M. (2001). Evolution of the insulin molecule: insights into structure-activity relationships. Peptides, 22(7), 1183-1193. PMID: 11445250.

