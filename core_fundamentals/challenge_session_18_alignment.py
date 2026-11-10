"""
Session 18 Challenge: Pairwise Sequence Alignment of HPV-16 and HPV-31 E7 Oncoproteins.
Author: Surena Ghadam
"""

from typing import Tuple
from Bio.Align import PairwiseAligner, substitution_matrices


def aligner_fix(
    mode: str = "global",
    substitution_matrix: str = "BLOSUM62",
    open_gap_score: float = -10.0,
    extend_gap_score: float = -0.5
) -> PairwiseAligner:
    """Configure and return a PairwiseAligner instance with specific scoring parameters."""
    if mode not in {"global", "local"}:
        raise ValueError(f"Invalid mode: {mode}. Must be 'global' or 'local'.")

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    aligner.substitution_matrix = substitution_matrices.load(substitution_matrix)
    
    return aligner


def calculate(seq1: str, seq2: str) -> Tuple[float, int, int]:
    """Calculate alignment metrics including identity percentage, length, and total gaps."""
    same_aa = 0
    gap_val = 0
    
    for char_q, char_t in zip(seq1, seq2):
        if char_q == char_t and char_q != "-":
            same_aa += 1
        elif char_q == "-" or char_t == "-":
            gap_val += 1

    length_aligned = len(seq1)
    identity = (same_aa / length_aligned) * 100.0 if length_aligned else 0.0
    
    return identity, length_aligned, gap_val


# Target E7 protein sequences (RefSeq: HPV-16 and HPV-31)
HPV_16 = (
    "MHGDTPTLHEYMLDLQPETTDLYCYEQLNDSSEEEDEIDGPAGQAEPDRAHYNIVTFCCKC"
    "DSTLRLCVQSTHVDIRTLEDLLMGTLGIVCPICSQKP"
)

HPV_31 = (
    "MRGETPTLQDYVLDLQPEATDLYCYEQLSDSSEEEDEIDGVDRPDGQAEQATSYKIVTCC"
    "CGCDSNVRLVVECTDGDIRQLQDLLLGTLNIVCPICAPKT"
)


def work() -> None:
    """Execute global and local alignments for HPV-16 and HPV-31 E7 proteins and report metrics."""
    # 1. Global Alignment (Needleman-Wunsch)
    aligner_global = aligner_fix(mode="global")
    alignment_g = aligner_global.align(HPV_16, HPV_31)
    global_best = alignment_g[0]
    score_g = global_best.score
    seq_q_g = global_best[0]
    seq_t_g = global_best[1]

    # 2. Local Alignment (Smith-Waterman)
    aligner_local = aligner_fix(mode="local")
    alignment_l = aligner_local.align(HPV_16, HPV_31)
    local_best = alignment_l[0]
    score_l = local_best.score
    seq_q_l = local_best[0]
    seq_t_l = local_best[1]

    # Calculate metrics
    identity_g, length_aligned_g, gap_val_g = calculate(seq1=seq_q_g, seq2=seq_t_g)
    identity_l, length_aligned_l, gap_val_l = calculate(seq1=seq_q_l, seq2=seq_t_l)

    # Output formatted results
    print(f"GLOBAL:     score: {score_g:<6.1f} | identity: {identity_g:5.2f}% | length_aligned: {length_aligned_g:<6} | gap Value: {gap_val_g:<6}")
    print(f"LOCAL:      score: {score_l:<6.1f} | identity: {identity_l:5.2f}% | length_aligned: {length_aligned_l:<6} | gap Value: {gap_val_l:<6}")


if __name__ == "__main__":
    work()
