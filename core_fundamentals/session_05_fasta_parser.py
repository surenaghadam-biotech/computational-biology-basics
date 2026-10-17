"""
Core Fundamental Module: FASTA Format Parser
Author: Surena Ghadam
"""
from pathlib import Path


def read_fasta(file_path: str) -> str:
    path = Path(file_path)
    seq = ""
    with open(path, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()
            if not line.startswith(">"):
                seq += line
    return seq


if __name__ == "__main__":
    file_name = "human_insulin.fasta"
    if Path(file_name).exists():
        insulin_seq = read_fasta(file_name)
        print(f"Loaded {len(insulin_seq)} bp")
