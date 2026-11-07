"""
"""
from Bio import Entrez, SeqIO
from utils import get_ncbi_email
from typing import List, Tuple
from Bio.SeqUtils import gc_fraction
import numpy as np
from matplotlib import pyplot as plt

#-------------------------------------------------------------------------------------------------
def window(
    seq_record:str, window_size:int, step_size:int=1)->List[Tuple[int,str,float]] :

    results: List[Tuple[int, str, float]] = []

    for start in range(0, len(seq_record) - window_size + 1, step_size):
        end = start + window_size
        window_kmer = seq_record[start:end]

        gc_percent = gc_fraction(window_kmer) * 100

        results.append((start, window_kmer, gc_percent))

    return results

#-------------------------------------------------------------------------------------------------
def plt_file(positions, gc_values, avr_gc, output_file):
    plt.figure(figsize=(10, 6))
    plt.plot(positions, gc_values, color='darkcyan', linewidth=1.5, marker='o', markersize=4)
    
    plt.title("GC Content per Window k-mer Sequence", fontsize=14, fontweight='bold')
    plt.ylabel("GC Content (%)", fontsize=12)
    plt.xlabel("Genomic Position (bp)", fontsize=12 )
    
    plt.ylim(0, 105)

    plt.grid(True, linestyle='--', alpha=0.5)

    plt.axhline(y=avr_gc, color='crimson', linestyle='--', label=f'Mean GC ({avr_gc:.1f}%)')
    plt.legend(loc='upper right')

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[+] Plot successfully saved to {output_file}")

#-------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    try:
        email = get_ncbi_email()
        Entrez.email = email
        with  Entrez.efetch(db="nuccore", id="NC_001526.4", rettype="fasta" , retmode="text") as handle:
            record = SeqIO.read(handle, "fasta")
    except Exception as e :
        print(f"network Error: {e}")
    else:

        seq_record = record.seq
        window_size = 200
        step_size = 50
        results = window(seq_record, window_size, step_size)
        positions = [r[0] for r in results]
        gc_values = [r[2] for r in results]
        avr_gc = np.mean(gc_values)

        output_file = r"D:\Computational_Biology\core_fundamentals\hpv16_gc_distribution.png"
        plt_file(positions,gc_values,avr_gc,output_file)

