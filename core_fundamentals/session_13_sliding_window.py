"""Session 13: Bioinformatics Algorithms - Sliding Window GC Content Analysis.

This module implements a sliding window algorithm to analyze local variations
in nucleotide composition (e.g., GC-content) across a genomic sequence.
"""

from typing import List, Tuple
from Bio.SeqUtils import gc_fraction # type: ignore[reportUnknownVariableType]


def calculate_sliding_window_gc(
    sequence: str, window_size: int, step_size: int = 1
) -> List[Tuple[int, str, float]]:
    """Scan a DNA sequence with a sliding window and compute GC percentage.

    Parameters:
        sequence (str): Target DNA sequence.
        window_size (int): Size of the sliding frame.
        step_size (int): Step step-forward step count (default: 1).

    Returns:
        List[Tuple[int, str, float]]: List of tuples containing:
            - Start index (0-based)
            - Window sub-sequence
            - GC content percentage (0.0 to 100.0)
    """

    results: List[Tuple[int, str, float]] = []
    seq_len: int = len(sequence)

    # Loop with a step from 0 until the window touches the end of sequence
    for start in range(0, seq_len - window_size + 1, step_size):
        end = start + window_size
        window_kmer = sequence[start:end]
        gc_score =float( gc_fraction(window_kmer) * 100.0)
        results.append((start, window_kmer, gc_score))

    return results


def visualize_terminal_sparkline(gc_results: List[Tuple[int, str, float]]) -> None:
    """Render a text-based visual bar chart in the terminal for GC content."""
    print("\n" + "=" * 60)
    print("      SLIDING WINDOW GC CONTENT PROFILER (ZJU LAB)")
    print("=" * 60)
    print(f"{'Pos':<6} | {'Subsequence':<12} | {'GC %':<6} | Visual Trend")
    print("-" * 60)

    for pos, kmer, gc in gc_results:
        # Create a text bar representing the percentage
        bar_length = int(gc // 5)  # Max length 20 characters
        bar = "█" * bar_length
        print(f"{pos:<6} | {kmer:<12} | {gc:>5.1f}% | {bar}")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    # Example: A sequence with an AT-rich region followed by a GC-rich region (like CpG island)
    mock_promoter_dna = "AAATATATATTAGCCGGCCGGCCATATAT"

    window: int = 6
    step: int = 2

    print(f"Analyzing sequence length: {len(mock_promoter_dna)} bp")
    print(f"Window Size: {window} bp | Step Size: {step} bp")

    # Run algorithm
    gc_profile = calculate_sliding_window_gc(
        mock_promoter_dna, window_size=window, step_size=step
    )

    # Visualize output
    visualize_terminal_sparkline(gc_profile)
