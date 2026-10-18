"""
Computational Biology Fundamentals - Session 11: In Silico Point Mutation Simulation
Author: Surena Ghadam
Context: Modeling the pathogenic Sickle Cell Disease mutation (HBB:c.20A>T, p.Glu6Val)
"""

from typing import Tuple
from Bio.Seq import Seq


def introduce_point_mutation(
    sequence: Seq,
    one_based_position: int,
    new_nucleotide: str
) -> Seq:
    """
    Introduces a single nucleotide variant (SNV) into a Seq object at a 1-based index.

    Parameters:
    -----------
    sequence : Seq
        The reference wild-type nucleotide sequence.
    one_based_position : int
        The 1-based genomic/cDNA coordinate of the mutation.
    new_nucleotide : str
        The replacement nucleotide base (A, T, C, or G).

    Returns:
    --------
    Seq
        The mutated sequence.
    """
    zero_based_index = one_based_position - 1
    seq_str = str(sequence)

    # Validate coordinate boundaries
    if not (0 <= zero_based_index < len(seq_str)):
        raise IndexError(f"Position {one_based_position} is out of bounds for sequence length {len(seq_str)}.")

    # Reconstruct the sequence with the mutated base
    mutated_str = (
        seq_str[:zero_based_index]
        + new_nucleotide.upper()
        + seq_str[zero_based_index + 1:]
    )
    return Seq(mutated_str)


def translate_cds_to_protein(cds_sequence: Seq) -> Seq:
    """
    Translates an in-frame CDS nucleotide sequence into an amino acid sequence.
    Stops at the first encountered stop codon.
    """
    # to_stop=True terminates translation at the first stop codon without raising warnings
    return cds_sequence.translate(to_stop=True)


def identify_amino_acid_variant(
    wt_protein: Seq,
    mut_protein: Seq
) -> Tuple[int, str, str]:
    """
    Compares wild-type and mutant proteins to locate the single amino acid substitution.

    Returns:
    --------
    Tuple[int, str, str]:
        (1-based position, WT residue, Mutant residue)
    """
    min_len = min(len(wt_protein), len(mut_protein))
    for i in range(min_len):
        if wt_protein[i] != mut_protein[i]:
            return i + 1, wt_protein[i], mut_protein[i]
    return 0, "", ""


def main() -> None:
    print("=" * 65)
    print("SESSION 11: IN SILICO HBB E6V MUTATION WORKSHOP")
    print("=" * 65)

    # 1. Authentic Human HBB CDS (Coding DNA Sequence: 444 nt -> 147 AA including Met)
    # Codon 1: ATG (M), Codon 2: GTG (V), Codon 3: CAC (H), Codon 4: CTG (L),
    # Codon 5: ACT (T), Codon 6: CCT (P), Codon 7: GAG (E) -> Sickle mutation site!
    # Note: Historically, HBB positions are often counted after initiator Met (p.Glu6Val)
    # In full CDS (HGVS nomenclature), this is c.20A>T (p.Glu7Val).
    
    # Let's take the first 30 nt of Human HBB CDS:
    # Pos: 1-3: ATG, 4-6: GTG, 7-9: CAC, 10-12: CTG, 13-15: ACT, 16-18: CCT, 19-21: GAG (Glu)
    wt_cds_prefix = Seq("ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTT")

    
    print(f"[*] Wild-Type CDS Prefix ({len(wt_cds_prefix)} nt): {wt_cds_prefix}")
    print(f"    WT Codon at pos 19-21       : {wt_cds_prefix[18:21]} -> Translates to: {wt_cds_prefix[18:21].translate()}")

    # 2. Simulate Point Mutation: Adenine at position 20 to Thymine (c.20A>T)
    # GAG -> GTG (Valine)
    mutated_cds_prefix = introduce_point_mutation(
        sequence=wt_cds_prefix,
        one_based_position=20,
        new_nucleotide="T"
    )

    print(f"[+] Mutated CDS Prefix (c.20A>T): {mutated_cds_prefix}")
    print(f"    Mutated Codon at pos 19-21  : {mutated_cds_prefix[18:21]} -> Translates to: {mutated_cds_prefix[18:21].translate()}")

    # 3. Translation comparison
    wt_prot = translate_cds_to_protein(wt_cds_prefix)
    mut_prot = translate_cds_to_protein(mutated_cds_prefix)

    print("\n" + "-" * 65)
    print("PROTEIN TRANSLATION COMPARISON")
    print("-" * 65)
    print(f"Wild-Type Peptide: {wt_prot}")
    print(f"Mutant Peptide   : {mut_prot}")

    # 4. Automated Variant Calling
    pos, wt_res, mut_res = identify_amino_acid_variant(wt_prot, mut_prot)
    if pos > 0:
        print(f"\n[!] Variant Detected: Position {pos} | {wt_res} -> {mut_res} (p.{wt_res}{pos}{mut_res})")
        print("    Biochemical Shift:")
        print(f"    - Wild-Type ({wt_res}): Glutamic Acid (Negatively charged, Hydrophilic)")
        print(f"    - Mutant    ({mut_res}): Valine (Non-polar, Hydrophobic -> triggers polymerization)")
    else:
        print("\n[-] No amino acid change detected (Synonymous mutation).")

    print("=" * 65)


if __name__ == "__main__":
    main()
