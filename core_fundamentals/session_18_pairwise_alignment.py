"""
Session 18: Pairwise Sequence Alignment & Substitution Matrices (BLOSUM62).

Comparative Evolutionary Analysis of Oncoprotein E7 across HPV-16 and HPV-18.
Demonstrates global and local sequence alignment using Biopython's modern
Bio.Align.PairwiseAligner API.

Author: Surena Ghadam (surenaghadam-biotech)
Environment: Python 3.8.10, Biopython 1.83
"""

from typing import Tuple
from Bio.Align import PairwiseAligner, substitution_matrices


def build_aligner(
    mode: str = "global",
    matrix_name: str = "BLOSUM62",
    open_gap_score: float = -10.0,
    extend_gap_score: float = -0.5,
) -> PairwiseAligner:
    """
    Configure and return a PairwiseAligner instance with affine gap penalties.

    :param mode: Alignment mode, either 'global' (Needleman-Wunsch)
                 or 'local' (Smith-Waterman).
    :param matrix_name: NCBI substitution matrix identifier (e.g., 'BLOSUM62').
    :param open_gap_score: Penalty for initiating a gap (affine model).
    :param extend_gap_score: Penalty for extending an existing gap.
    :return: Fully configured PairwiseAligner object.
    :raises ValueError: If an unsupported mode is provided.
    """
    if mode not in {"global", "local"}:
        raise ValueError(f"Invalid mode '{mode}'. Expected 'global' or 'local'.")

    aligner = PairwiseAligner()
    aligner.mode = mode

    # Load amino acid substitution matrix from Biopython repository
    matrix = substitution_matrices.load(matrix_name)
    aligner.substitution_matrix = matrix

    # Set affine gap penalties (scores are negative values)
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    return aligner


def calculate_identity_metrics(
    aligned_query: str, aligned_target: str
) -> Tuple[int, int, float]:
    """
    Calculate exact matches, total aligned positions, and sequence identity percentage.

    :param aligned_query: The query sequence with gap characters ('-').
    :param aligned_target: The target sequence with gap characters ('-').
    :return: Tuple containing (identical_positions, total_columns, identity_percentage).
    """
    identical_positions = 0
    total_columns = len(aligned_query)

    for char_q, char_t in zip(aligned_query, aligned_target):
        if char_q == char_t and char_q != "-":
            identical_positions += 1

    identity_pct = (
        (identical_positions / total_columns) * 100.0 if total_columns > 0 else 0.0
    )
    return identical_positions, total_columns, identity_pct


def run_comparative_analysis() -> None:
    """Execute pairwise alignment between HPV-16 E7 and HPV-18 E7 oncoproteins."""
    # HPV-16 E7 oncoprotein (Reference: NC_001526.4 / P03129, 98 aa)
    hpv16_e7: str = (
        "MHGDTPTLHEYMLDLQPETTDLYCYEQLNDSSEEEDEIDGPAGQAEPDRAHYNIVTFCCKC"
        "DSTLRLCVQSTHVDIRTLEDLLMGTLGIVCPICSQKP"
    )

    # HPV-18 E7 oncoprotein (Reference: NC_001357.1 / P06788, 105 aa)
    hpv18_e7: str = (
        "MHGPKATLQDIVLHLEPQNEIPVDLLCHEQLSDSEEENDEIDGVNHQHLPARRAEPQRHTML"
        "CMCCKCEARIELVVESSADDLRAFQQLFLNTLSFVCPWCASQQ"
    )

    print("=" * 70)
    print("      HPV-16 E7 vs HPV-18 E7 PAIRWISE SEQUENCE ALIGNMENT (BLOSUM62)")
    print("=" * 70)

    # 1. Initialize Aligner in Global Mode (Needleman-Wunsch)
    aligner = build_aligner(
        mode="global",
        matrix_name="BLOSUM62",
        open_gap_score=-10.0,
        extend_gap_score=-0.5,
    )

    # 2. Compute alignment score and get optimal alignments
    score: float = aligner.score(hpv16_e7, hpv18_e7)
    alignments = aligner.align(hpv16_e7, hpv18_e7)

    # Take the top optimal alignment
    best_alignment = alignments[0]

    # Convert alignment lines to strings for metric extraction
    # best_alignment[0] is query (with gaps), best_alignment[1] is target (with gaps)
    aligned_seq1 = best_alignment[0]
    aligned_seq2 = best_alignment[1]

    matches, cols, identity = calculate_identity_metrics(aligned_seq1, aligned_seq2)

    # 3. Print Clinical & Evolutionary Metrics
    print(f"Algorithm Mode         : Global (Needleman-Wunsch)")
    print(f"Substitution Matrix    : BLOSUM62")
    print(f"Affine Gap Penalties   : Open=-10.0, Extend=-0.5")
    print(f"Optimal Alignment Score: {score:.2f}")
    print(f"Total Aligned Columns  : {cols}")
    print(f"Identical Residues     : {matches} / {cols}")
    print(f"Sequence Identity      : {identity:.2f}%")
    print("-" * 70)
    print("Alignment Visual Trace:")
    print(best_alignment)
    print("=" * 70)


if __name__ == "__main__":
    run_comparative_analysis()
