import sys

print("=" * 60)
print(f"Python Version: {sys.version.split()[0]}")
print("=" * 60)

# تست بیوپایتون به شکل استاندارد
try:
    import Bio
    from Bio.Seq import Seq
    print(f"[OK] Biopython          | نسخه: {Bio.__version__} | پکیج مرجع بیوانفورماتیک")
    
    # تست رونویسی و ترجمه ژنتیکی در حافظه
    sample_dna = Seq("ATGCGATCGATCGATCGATAG")
    mrna = sample_dna.transcribe()
    protein = sample_dna.translate()
    
    print("-" * 60)
    print("تست محاسباتی ژنتیک مولکولی (In Silico Validation):")
    print(f"   - توالی DNA نمونه:      {sample_dna}")
    print(f"   - رونویسی (Transcription): {mrna}")
    print(f"   - ترجمه (Translation):     {protein}")
    print("=" * 60)
    print(" وضعیت: تبریک! سیستم ۱۰۰٪ آماده ورود به جلسه اول (فاز ۱) است.")

except ImportError as e:
    print(f"[ERROR] Biopython هنوز نصب نشده: {e}")
import sys

print("=" * 60)
print(f"Python Version: {sys.version.split()[0]}")
print("=" * 60)

# تست بیوپایتون به شکل استاندارد
try:
    import Bio
    from Bio.Seq import Seq
    print(f"[OK] Biopython          | نسخه: {Bio.__version__} | پکیج مرجع بیوانفورماتیک")
    
    # تست رونویسی و ترجمه ژنتیکی در حافظه
    sample_dna = Seq("ATGCGATCGATCGATCGATAG")
    mrna = sample_dna.transcribe()
    protein = sample_dna.translate()
    
    print("-" * 60)
    print("تست محاسباتی ژنتیک مولکولی (In Silico Validation):")
    print(f"   - توالی DNA نمونه:      {sample_dna}")
    print(f"   - رونویسی (Transcription): {mrna}")
    print(f"   - ترجمه (Translation):     {protein}")
    print("=" * 60)
    print(" وضعیت: تبریک! سیستم ۱۰۰٪ آماده ورود به جلسه اول (فاز ۱) است.")

except ImportError as e:
    print(f"[ERROR] Biopython هنوز نصب نشده: {e}")
