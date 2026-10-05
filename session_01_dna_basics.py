"""
پروژه بیوانفورماتیک محاسباتی - جلسه اول
موضوع: تحلیل پایه‌ای توالی DNA، اسلایسینگ و محاسبه GC-Content
نویسنده: سورنا قدم
"""

# ۱. تعریف متغیرها (توالی ژن هدف)
# دقت کن که رشته‌ها داخل کوتیشن ' ' قرار می‌گیرند.
gene_name = "LONGEVITY_CANDIDATE_GENE_1"
dna_seq = "ATGCGATCGATCGCTAGCTAGCTAGCGCTAGCTAGCTAATCGATCGATCGATAG"

print("=" * 65)
print(f"🧬 آنالیز بیوانفورماتیک برای ژن: {gene_name}")
print("=" * 65)

# ۲. به دست آوردن طول توالی با استفاده از تابع پیش‌فرض len()
seq_length = len(dna_seq)
print(f"1. طول توالی ژن: {seq_length} جفت‌باز (bp)")

# ۳. دسترسی مستقیم به نوکلئوتیدها (Indexing)
first_base = dna_seq[0]  # اولین نوکلئوتید (ایندکس ۰)
print(f"2. اولین باز در سمت 5-پرایم ('5): {first_base}")

# ۴. استخراج زیرتوالی یا برش (Slicing) با ساختار: [شروع:پایان]
# ایندکس پایان خوانده نمی‌شود (Exclusive است)
start_codon = dna_seq[0:3]  # کاراکترهای شماره ۰، ۱ و ۲
print(f"3. کدون آغازین استخراج‌شده: {start_codon}")

# ۵. خواندن از انتهای توالی با اندیس منفی (3 کاراکتر آخر)
stop_codon_candidate = dna_seq[-3:]
print(f"4. کاندید کدون پایانی در سمت 3-پرایم ('3): {stop_codon_candidate}")

# ۶. شمارش فراوانی نوکلئوتیدها با متد .count()
count_a = dna_seq.count("A")
count_t = dna_seq.count("T")
count_c = dna_seq.count("C")
count_g = dna_seq.count("G")

print("-" * 65)
print("5. فراوانی بازهای نیتروژنی:")
print(f"   - آدنین (A): {count_a} عدد")
print(f"   - تیمین (T): {count_t} عدد")
print(f"   - سیتوزین (C): {count_c} عدد")
print(f"   - گوانین (G): {count_g} عدد")

# ۷. محاسبه بیولوژیکی GC-Content (فرمول: مجموع G و C تقسیم بر کل طول ضرب در ۱۰۰)
gc_count = count_g + count_c
gc_content = (gc_count / seq_length) * 100

print("-" * 65)
# با استفاده از 2f. عدد اعشاری را تا دو رقم گرد می‌کنیم
print(f"6. درصد GC-Content توالی: {gc_content:.2f}%")

if gc_content > 50.0:
    print("   [تحلیل زیستی: توالی غنی از GC است؛ پایداری ساختاری و Tm بالا]")
else:
    print("   [تحلیل زیستی: توالی غنی از AT است؛ احتمال دناتوراسیون در دمای پایین‌تر]")

print("=" * 65)
