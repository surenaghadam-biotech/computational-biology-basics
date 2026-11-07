"""
Session 15: Matplotlib-based GC Content Visualizer.
Author: Surena Ghadam
Description: Generates publication-ready genomic GC content distribution plots.
"""

from typing import List, Tuple
import matplotlib.pyplot as plt

from session_13_sliding_window import calculate_sliding_window_gc


def plot_gc_content(results: List[Tuple[int, str, float]], output_file: str) -> None:
    """
    Plots the GC content distribution and saves as a PNG file.
    
    Args:
        results: List of (position, kmer, gc_percentage).
        output_file: Path to save the plot.
    """
    # positions = [r[0] for r in results]
    positions = [r[0] for r in results]
    gc_values = [r[2] for r in results]

    plt.figure(figsize=(10, 6))
    plt.plot(positions, gc_values, color='teal', linewidth=1.5, marker='o', markersize=4)
    plt.title("GC Content Distribution across Genomic Sequence", fontsize=14, fontweight='bold')
    plt.xlabel("Position (bp)", fontsize=12)
    plt.ylabel("GC Content (%)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # تنظیم محدودیت محور Y برای خوانایی بهتر (از 0 تا 100 درصد)
    plt.ylim(0, 105)
    
    # plt.xticks(positions) 

    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[+] Plot successfully saved to {output_file}")


if __name__ == "__main__":
    # توالی تستی نمونه
    test_sequence = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    window_size = 6
    step_size = 2
    
    # محاسبه پنجره لغزان
    gc_results = calculate_sliding_window_gc(test_sequence, window_size, step_size)
    
    # مسیر خروجی امن با Raw String و پسوند png
    output_path = r"D:\Computational_Biology\core_fundamentals\gc_content_plot.png"
    
    # رسم و ذخیره نمودار
    plot_gc_content(gc_results, output_path)
