
from typing import List, Tuple
import matplotlib.pyplot as plt

from session_13_sliding_window import calculate_sliding_window_gc

def plot_gc_content(results: List[Tuple[int, str, float]], output_file: str) -> None:

    # ۱. استخراج داده‌ها به شیوه درست
    positions = [r[0] for r in results]      # مقادیر عددی برای حفظ ساختار پیوسته
    sequences = [r[1] for r in results]      # متن توالی‌ها برای درج روی خط‌کش
    gc_values = [r[2] for r in results]      # مقادیر درصد GC

    # ۲. رسم نمودار با محور عددی
    plt.figure(figsize=(12, 6))
    plt.plot(positions, gc_values, color='teal', linewidth=1.5, marker='o', markersize=4)

    # ۳. جادوی کار: نگاشت متن توالی‌ها روی نقاط دقیق خط‌کش
    # پارامتر اول: مکان اعداد (ticks) | پارامتر دوم: برچسب‌هایی که چاپ می‌شوند (labels)
    plt.xticks(ticks=positions, labels=sequences, rotation=45, ha='right', fontsize=9)

    # ۴. حاشیه‌ها و نمایش
    plt.title("GC Content per Window k-mer Sequence", fontsize=14, fontweight='bold')
    plt.xlabel("Window k-mer Sequence", fontsize=12)
    plt.ylabel("GC Content (%)", fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.ylim(0, 105)
    plt.tight_layout()  # جلوگیری از بریده شدن حروف کج‌شده در لبه تصویر


    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[+] Plot successfully saved to {output_file}")


if __name__ == "__main__":
    # توالی تستی نمونه
    test_sequence = "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC"
    window_size = 3
    step_size = 3
    
    # محاسبه پنجره لغزان
    gc_results = calculate_sliding_window_gc(test_sequence, window_size, step_size)
    
    # مسیر خروجی امن با Raw String و پسوند png
    output_path = r"D:\Computational_Biology\core_fundamentals\gc_content_plot_str1.png"
    
    # رسم و ذخیره نمودار
    plot_gc_content(gc_results, output_path)
