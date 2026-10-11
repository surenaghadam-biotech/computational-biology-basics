def calculate_identity_score(seq1, seq2):
    """
    Calculates the sequence identity percentage between two sequences
    based on the overlapping length (min_len).
    """
    min_len = min(len(seq1), len(seq2))
    same_position = 0
    
    for i in range(min_len):
        if seq1[i] == seq2[i]:
            same_position += 1
            
    identity = (same_position / min_len) * 100
    return identity
