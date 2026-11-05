"""
Session 12: Sequence Alignment Fundamentals - Dot Plot Analysis
Author: Surena Ghadam
Date: 2026

This script implements a 2D Dot Plot Matrix from scratch to visualize
sequence homology and local alignments between two nucleotide sequences.
"""

import numpy as np


def build_dotplot_matrix(seq1: str, seq2: str) -> np.ndarray:
    """
    Constructs a binary dot plot matrix comparing two biological sequences.
    
    Args:
        seq1: First sequence string (X-axis).
        seq2: Second sequence string (Y-axis).
        
    Returns:
        np.ndarray: 2D binary numpy array where 1 represents a match.
    """
    matrix = np.zeros((len(seq2), len(seq1)), dtype=int)
    
    for row_idx, char2 in enumerate(seq2):
        for col_idx, char1 in enumerate(seq1):
            if char1 == char2:
                matrix[row_idx, col_idx] = 1
                
    return matrix


def display_ascii_dotplot(seq1: str, seq2: str, matrix: np.ndarray) -> None:
    """
    Prints a clean, formatted ASCII dot plot to the terminal.
    
    Args:
        seq1: First sequence string (columns).
        seq2: Second sequence string (rows).
        matrix: 2D binary numpy matrix.
    """
    print("\n" + "=" * 50)
    print("       BIOINFORMATICS TERMINAL DOT PLOT")
    print("=" * 50)
    
    # Print header with columns
    header = "    " + " ".join(list(seq1))
    print(header)
    print("   +" + "-" * (len(seq1) * 2))

    # Print each row with its corresponding character
    for row_idx, char2 in enumerate(seq2):
        row_symbols = []
        for col_idx in range(len(seq1)):
            if matrix[row_idx, col_idx] == 1:
                row_symbols.append("●")  # Match symbol
            else:
                row_symbols.append("·")  # Mismatch symbol
        print(f" {char2} | " + " ".join(row_symbols))
        
    print("=" * 50 + "\n")


def calculate_identity_percentage(matrix: np.ndarray, seq1: str, seq2: str) -> float:
    """
    Calculates overall identity percentage along the main diagonal.
    
    Args:
        matrix: 2D numpy dot plot matrix.
        seq1: First sequence string.
        seq2: Second sequence string.
        
    Returns:
        float: Percentage of identical bases.
    """
    min_len = min(len(seq1), len(seq2))
    diagonal_matches = np.trace(matrix[:min_len, :min_len])
    identity = (diagonal_matches / min_len) * 100
    return float(identity)


def main() -> None:
    # Example 1: Comparing wild-type vs point mutated DNA fragment
    seq_wildtype = "ATGGTGCACCTGACTCCT"
    seq_mutated  = "ATGGTGCATCTGACTCCT"  # C -> T mutation at position 9
    
    print(f"Seq 1 (Wildtype): {seq_wildtype}")
    print(f"Seq 2 (Mutated) : {seq_mutated}")
    
    # 1. Build matrix
    dot_matrix = build_dotplot_matrix(seq_wildtype, seq_mutated)
    
    # 2. Print ASCII visual representation
    display_ascii_dotplot(seq_wildtype, seq_mutated, dot_matrix)
    
    # 3. Calculate sequence identity
    identity_score = calculate_identity_percentage(dot_matrix, seq_wildtype, seq_mutated)
    print(f"Main Diagonal Sequence Identity: {identity_score:.2f}%")


if __name__ == "__main__":
    main()
