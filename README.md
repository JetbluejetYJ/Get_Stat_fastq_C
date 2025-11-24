# 📊 Get_stat_fastq with C (kseq + zlib + pthreads)

A lightweight, multi-threaded C utility that scans gzipped FASTQ files and emits compact QC summaries per sample. It uses `kseq.h` for FASTQ parsing, `zlib` for streaming decompression, and POSIX threads to process multiple sample prefixes concurrently.

---

## ✨ Features

- 🔎 Discovers samples by scanning a directory for files matching `*_1.fastq.gz`.
- 🔗 Treats `<prefix>_1.fastq.gz` and `<prefix>_2.fastq.gz` as a paired-end set; if R2 is absent it processes R1 as single-end.
- 🚀 One thread per sample prefix (up to 128 concurrently).
- 📊 Per-sample output (`.sqs`) with:
  - Total bases and total reads
  - Average read size
  - Base composition counts (A, C, G, T, N)
  - GC percentage, N percentage
  - Q20/Q30 base counts and percentages
  - Per-mate (R1, R2) breakdown when paired-end is present
- 💾 Large zlib buffer (8 MB) via `gzbuffer` to improve I/O throughput.
- 🔁 Skip-if-present behavior: if `<prefix>.sqs` already exists and is non-empty, that sample is skipped.

---

## 🧭 How it works

1. Directory scan
   - The program expects a single argument: the directory containing gzipped FASTQs.
   - It looks for filenames ending with `_1.fastq.gz`. The substring before `_1.fastq.gz` is treated as the sample prefix.
   - If `<prefix>_2.fastq.gz` exists, the sample is considered paired-end.

2. Concurrency
   - Each discovered prefix is processed in its own thread (maximum 128). Threads are joined at the end.

3. FASTQ parsing
   - `kseq.h` is used on top of `gzFile` (`gzopen`, `gzread`) to stream reads.
   - For each read, the sequence and quality arrays are traversed once.

4. Output
   - Results are written to `<dir>/<prefix>.sqs`.
   - A top summary row is followed by detailed sections. If paired, R1 and R2 sections are emitted.

---

## 🧮 Calculations and formulas

Let:
- $B$ = total number of bases across all processed reads.
- $R$ = total number of reads processed.
- $L_r$ = length of read $r$.
- $A, C, G, T, N$ = base counts (uppercase only; other letters are not tallied as $N$ in this implementation).
- $B_{Q\ge20}$, $B_{Q\ge30}$ = counts of bases with Phred quality $Q$ meeting the thresholds.

Then:
- Total bases: $B = \sum_{r=1}^{R} L_r$
- Average read size: $\overline{L} = \dfrac{\sum_{r=1}^{R} L_r}{\max(R,1)}$
- GC percentage: $\mathrm{GC\%} = 100 \times \dfrac{G + C}{B}$
- N percentage: $\mathrm{N\%} = 100 \times \dfrac{N}{B}$
- Q20 percentage: $\mathrm{Q20\%} = 100 \times \dfrac{B_{Q\ge 20}}{B}$
- Q30 percentage: $\mathrm{Q30\%} = 100 \times \dfrac{B_{Q\ge 30}}{B}$



Quality handling:
- Qualities are treated as Sanger/Illumina Phred+33. The code compares the ASCII code of each quality character to the threshold offset: it increments $B_{Q\ge20}$ if `qv >= 20 + 33`, and $B_{Q\ge30}$ if `qv >= 30 + 33`. In other words, $Q \ge 20$ or $Q \ge 30$ inclusively.
- No recalibration or trimming is performed.

Per-mate (paired-end):
- The program maintains separate counters for R1 and R2: $B^{(1)}$, $B^{(2)}$, $R^{(1)}$, $R^{(2)}$, base counts, and Q20/Q30 counts. The same formulas apply to each mate using its own denominators.

Notes and caveats:
- Base counting considers only uppercase 'A', 'C', 'G', 'T', 'N'. Lowercase bases or other IUPAC letters are not counted (and are not redirected to $N$).
- When paired files differ in the number of records, processing stops when either file reaches EOF; any trailing reads in the longer file are not processed.
- The first summary line percentages are rounded via simple floating-point formatting in the code.

---

## 📄 Output format

For each sample, the `.sqs` file contains:

1) A tab-separated summary line:
```
<SampleName>    <TotalBase>    <TotalRead>    <N%>    <GC%>    <Q20%>    <Q30%>
```

2) Named blocks with detailed counts (overall first). If paired, additional R1 and R2 blocks follow.
Example (illustrative numbers):
```
SAMPLE_X    150000000    1000000    0.1200    41.50    97.80    92.25
SampleName : SAMPLE_X
Total A : 45000000
Total C : 30000000
Total G : 32250000
Total T : 42250000
Total N : 180000
Q30 Bases : 138375000
Q20 Bases : 146700000
Avg.ReadSize : 150.00

SAMPLE_X_R1    75000000    500000    0.1100    41.60    97.90    92.40
...
SAMPLE_X_R2    75000000    500000    0.1300    41.40    97.70    92.10
...
```

---

## 🛠️ Build

Dependencies:
- `zlib` (for `gzopen/gzread`)
- `pthread` (POSIX threads)
- `kseq.h` (single-header FASTA/FASTQ parser by Heng Li). Place `kseq.h` next to your source or add an include directory with `-I`.

Example :
```bash
gcc -O3 -std=c11 -Wall -Wextra -pthread -o get_stat_fastq_C get_stat_fastq_C.c -lz
```

If `kseq.h` is not in the same directory:
```bash
gcc -O3 -std=c11 -Wall -Wextra -pthread -I/path/to/klib -o get_stat_fastq_C get_stat_fastq_C.c -lz
```

Performance tips:
- Consider `-march=native` on bare-metal environments.
- Ensure your storage is not the bottleneck; this tool streams gzip data and benefits from SSDs/NVMe.

---

## 🚀 Usage

```text
Usage: <get_stat_fastq_C> [path of order data]
```

Examples:
```bash
# Process all samples under a directory
./get_stat_fastq_C /data/runs/2025-11-01/

# Skip already summarized samples (existing non-empty .sqs files)
./get_stat_fastq_C /data/runs/2025-11-01/
```

Behavioral details:
- Samples are discovered by `_1.fastq.gz` naming.
- If `<prefix>.sqs` exists and is non-empty, that sample is skipped.
- Up to 128 samples per run (thread slots are fixed at 128 in the code).

---

## ⚠️ Limitations

- Only gzipped FASTQ files are supported.
- Maximum concurrent samples: 128 (as coded).
- Quality scale is assumed to be Phred+33; Phred+64 data will be misinterpreted.
- Base tallies ignore lowercase and non-ACGT IUPAC codes (not mapped to `N`).
- When R1/R2 lengths differ (record counts), the loop stops at the shorter file’s EOF.

---

## 🧪 Validation checklist

- Confirm your FASTQs use Phred+33 (modern Illumina).
- Check that sample names follow `<prefix>_1.fastq.gz` and `<prefix>_2.fastq.gz`.
- Verify that `.sqs` is created in the same directory as the inputs.

---

## 📚 Acknowledgments

- `kseq.h` by Heng Li.
- `zlib` for gzip streaming.

