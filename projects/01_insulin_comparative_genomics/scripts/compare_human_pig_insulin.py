"""
Research Project: Comparative Insulin Sequence Analysis (Human vs. Pig)
PubMed Reference: PMID: 11445250
Author: Surena Ghadam
"""
from pathlib import Path
from core_fundamentals.session_05_fasta_parser import read_fasta
from core_fundamentals.def_dna_to_protein import translate_dna
from calculate_identity_score import calculate_identity_score
from def_dotplot_alignment import generate_dot_plot

PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_DIR / "data"

pig_mrna = read_fasta(str(DATA_DIR / "pig_insulin.fasta"))
pig_pro = translate_dna(pig_mrna)

human_mrna = read_fasta(str(DATA_DIR / "human_insulin.fasta"))
human_pro = translate_dna(human_mrna)

generate_dot_plot(pig_pro, human_pro)
score = calculate_identity_score(pig_pro, human_pro)

print(f"Identity Score: {score:.2f}%")
