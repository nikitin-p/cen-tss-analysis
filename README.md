# cen-tss-analysis

Analysis of transcription start sites in human alpha-satellite pericentromeric regions.

This pipeline takes CAGE BAM files that were mapped to the T2T reference genome with STAR,
and produces strand-specific 5'-end coverage bigWig files suitable for TSS analysis and
initiator dinucleotide extraction.

---

## Pipeline overview

The pipeline has six steps. Each step is a standalone script in `scripts/` and can be
run independently if needed.

```text
STAR BAM (input)
    │
    ▼ scripts/get_remap.sh
Reads from centromeric regions + unmapped reads → merged FASTQ
    │
    ▼ scripts/bowtie2_remap.sh
Bowtie2 BAM (multimappers retained, local alignment)
    │
    ▼ scripts/select_tss_alignments.sh  (uses scripts/pyselectal.py)
BAM filtered to reads with exactly one 5'-soft-clipped G
    │
    ▼ scripts/xm_filter.sh
BAM filtered to reads with < 10% mismatches (XM tag)
    │
    ▼ scripts/filter_transposons.sh
BAM with transposon-overlapping reads removed
    │
    ▼ scripts/bamcoverage_5p.sh
SAMPLE.fwd.bw  (forward strand, + )
SAMPLE.rev.bw  (reverse strand, - )
```

---

## Repository structure

```text
cen-tss-analysis/
├── run_pipeline.sh               Main orchestration script
├── pipeline.conf                 Configuration: paths, parameters, thread count
├── scripts/
│   ├── get_remap.sh              Step 1: extract centromeric and unmapped reads
│   ├── bowtie2_remap.sh          Step 2: remap with bowtie2
│   ├── select_tss_alignments.sh  Step 3: filter to 5'-G soft-clipped reads
│   ├── pyselectal.py             Python tool used by select_tss_alignments.sh
│   ├── xm_filter.sh              Step 4: filter by mismatch rate
│   ├── filter_transposons.sh     Step 5: remove transposon-overlapping reads
│   └── bamcoverage_5p.sh         Step 6: generate strand-specific bigWigs
├── scripts_old/                  Original scripts kept for reference (not used)
└── README.md
```

### What each file does

| File | Purpose |
| ---- | ------- |
| `run_pipeline.sh` | Orchestrates all six steps, handles cleanup, logging, resume, and dry-run |
| `pipeline.conf` | All paths and tunable parameters; edit this before running |
| `scripts/get_remap.sh` | Extracts reads overlapping centromeric BED regions and all unmapped reads from the input STAR BAM, merges them into a FASTQ for remapping |
| `scripts/bowtie2_remap.sh` | Remaps the FASTQ to the centromeric reference with bowtie2 in local mode, retaining multimappers |
| `scripts/select_tss_alignments.sh` | Keeps only reads with exactly one 5'-soft-clipped base that is G (forward) or C (reverse); delegates to pyselectal.py |
| `scripts/pyselectal.py` | Python BAM filter for 5'-end soft-clip and sequence selection; checks both CIGAR pattern and actual nucleotide identity |
| `scripts/xm_filter.sh` | Keeps reads where (XM tag / read length) < threshold (default < 10% mismatches) |
| `scripts/filter_transposons.sh` | Removes any alignment that overlaps the centromeric transposon BED using `bedtools intersect -v` |
| `scripts/bamcoverage_5p.sh` | Runs bamCoverage twice (once per strand) to produce 1 bp resolution 5'-end bigWig files |
| `scripts_old/` | Original scripts from before the refactor; kept as reference but not executed by the pipeline |

---

## Testing

A set of synthetic test BAMs and BEDs is provided to verify steps 1, 3, 4, 5, and 6
without needing real CAGE data or a bowtie2 index.

### Generate test data

```bash
python tests/make_test_data.py
```

Requires `pysam` and `samtools`. Writes the following to `tests/data/`:

| File | Purpose |
| ---- | ------- |
| `centromere.bed` | `chr1:500-2000` — centromere regions for step 1 |
| `transposon.bed` | `HOR_chr1:1000-2000` — transposon regions for step 5 |
| `step1_input.bam` | 5 centromeric + 3 unmapped + 4 off-target reads (STAR-like) |
| `step3_input.bam` | 3 reads that pass 5'-G filter + 4 that fail |
| `step4_input.bam` | 2 pass XM filter, 2 fail, 1 dropped (no XM tag) |
| `step5_input.bam` | 2 outside transposon (pass) + 2 overlapping (fail) |

### Run tests

```bash
bash tests/run_tests.sh
```

Requires `samtools`, `bedtools`, `bamCoverage` (deepTools). Steps 2 (bowtie2 remap)
is not tested as it requires a bowtie2 index.

Expected output:

```text
[PASS] reads in FASTQ (cen+unmapped)              expected 8, got 8
[PASS] reads passing 5'-G filter                  expected 3, got 3
[PASS] reads passing XM filter                    expected 2, got 2
[PASS] reads passing transposon filter            expected 2, got 2
[PASS] step6_output.fwd.bw                        exists and non-empty
[PASS] step6_output.rev.bw                        exists and non-empty
────────────────────────────────────────────────────────
Results: 6 / 6 passed
```

---

## Dependencies

| Tool | Version | Notes |
| ---- | ------- | ----- |
| bash | ≥ 4.2 | Array syntax required |
| samtools | ≥ 1.13 | |
| bedtools | ≥ 2.29 | |
| bowtie2 | ≥ 2.4 | |
| bamCoverage | ≥ 3.5 (deepTools) | |
| python3 | ≥ 3.8 | |
| pysam | ≥ 0.19 | Python package; used by pyselectal.py |

Install deepTools and pysam via conda or pip:

```bash
conda install -c bioconda deeptools pysam
# or
pip install deeptools pysam
```

---

## Configuration

Edit `pipeline.conf` before running. The variables are:

| Variable | Description |
| -------- | ----------- |
| `CENTROMERE_BED` | BED file of centromeric regions (HOR/alpha-satellite) |
| `TRANSPOSON_BED` | BED file of centromeric transposable elements |
| `BOWTIE2_INDEX_PREFIX` | Prefix of the bowtie2 index for the centromeric reference |
| `BOWTIE2_EXTRA_ARGS` | Bash array of bowtie2 flags (default: `-k 800 -D 1600 -R 1 --local`) |
| `PYSELECTAL_ARGS` | Bash array of pyselectal.py flags (default: `-n 1 -m 1 -x G`) |
| `PYSELECTAL_PY` | Path to `pyselectal.py`; defaults to `scripts/pyselectal.py` in this repo |
| `MIN_MATCH_FRAC` | Minimum match fraction for XM filter (default: `0.90` = < 10% mismatches) |
| `THREADS` | Default thread count; overridden at runtime with `-t` |

---

## Usage

### Single BAM file

```bash
bash run_pipeline.sh \
  -i data/sampleA.bam \
  -o results/ \
  -c pipeline.conf
```

### Directory of BAM files

```bash
bash run_pipeline.sh \
  -i data/bams/ \
  -o results/ \
  -c pipeline.conf \
  -t 16
```

All `.bam` files directly inside `data/bams/` are processed. Subdirectories are not scanned.

### Resume an interrupted run

```bash
bash run_pipeline.sh -i data/bams/ -o results/ --resume
```

Steps whose output files already exist are skipped. Useful after a crash or when adding
new samples to an existing output directory.

### Keep intermediate files

```bash
bash run_pipeline.sh -i data/bams/ -o results/ --keep-intermediates
```

By default all intermediate BAM and FASTQ files are deleted immediately after the next
step completes successfully (see Cleanup policy below).

### Dry run

```bash
bash run_pipeline.sh -i data/bams/ -o results/ --dry-run
```

Prints all commands that would be executed without running anything. Config and tool
validation still runs, so you can catch path errors before submitting a long job.

---

## Output files

For each input BAM named `SAMPLE.bam`, the pipeline creates `results/SAMPLE/`:

| File | Description |
| ---- | ----------- |
| `SAMPLE.fwd.bw` | 5'-end coverage bigWig, forward strand (+) |
| `SAMPLE.rev.bw` | 5'-end coverage bigWig, reverse strand (−) |
| `SAMPLE.log` | Full log of all tool output for this sample |

If `--keep-intermediates` is passed, these files are also retained:

| File | Step that creates it |
| ---- | -------------------- |
| `SAMPLE.remap.fastq.gz` | Step 1 (get_remap) |
| `SAMPLE.bowtie2.bam` | Step 2 (bowtie2_remap) |
| `SAMPLE.softclip_g.bam` | Step 3 (select_tss_alignments) |
| `SAMPLE.xm_filtered.bam` | Step 4 (xm_filter) |
| `SAMPLE.no_transposons.bam` | Step 5 (filter_transposons) |

---

## Sample naming

The sample name is derived by stripping the `.bam` extension from the input filename:

```text
/path/to/sampleA.bam          → sampleA
/path/to/A12_rep2.sorted.bam  → A12_rep2.sorted
```

No assumptions are made about internal naming conventions (no underscore splitting,
no suffix removal beyond `.bam`). Output directories and filenames use this name verbatim.

---

## Cleanup policy

Intermediate files are deleted **immediately after the downstream step succeeds** and
its output has been validated as non-empty. The deletion sequence is:

1. After step 2 completes: delete `SAMPLE.remap.fastq.gz`
2. After step 3 completes: delete `SAMPLE.bowtie2.bam` and its index
3. After step 4 completes: delete `SAMPLE.softclip_g.bam` and its index
4. After step 5 completes: delete `SAMPLE.xm_filtered.bam` and its index
5. After step 6 completes: delete `SAMPLE.no_transposons.bam` and its index

`trap`-based cleanup is intentionally **not** used for intermediates, so that a failed
run leaves the last successfully created file on disk for debugging.

Use `--keep-intermediates` to disable all deletion.

---

## Assumptions and known ambiguities inherited from the original scripts

### 1. Input to bowtie2 (step 2)

Both centromeric reads (mapped to the centromeric BED regions) **and** unmapped reads
are extracted from the STAR BAM and combined into one FASTQ. Both are remapped together.
This matches the intent of the original pipeline design.

### 2. 5'-G base check (step 3) — behaviour change from original

The original `get5GsoftclipBAM.sh` selected reads by CIGAR pattern only (`1S` prefix),
without checking whether the soft-clipped base was actually G. `select_tss_alignments.sh`
uses `pyselectal.py -n 1 -m 1 -x G`, which checks **both** the CIGAR and the nucleotide
identity. This is biologically more precise for CAGE TSS identification and is consistent
with the intent indicated by the `pyselectal.py` call in the original pipeline description.
Reads where the 5'-soft-clipped base is not G will be excluded here but were retained
by the old script.

### 3. XM mismatch threshold

The original `xm_filter.sh` used `XM/read_length < 0.10` (strict less-than).
The refactored version uses `<= (1 - MIN_MATCH_FRAC)`, which with the default `0.90`
gives `<= 0.10` (less-than-or-equal). The difference affects only reads with exactly
10% mismatches, which is negligible in practice.

### 4. Transposon filter

Any alignment with **any** overlap to the transposon BED is removed. Partial overlaps
are treated the same as full overlaps. The original intent was to remove transposon-
derived signal, and `bedtools intersect -v` with default parameters implements this.

### 5. BigWig position and the 1 bp soft-clip offset

Coverage is reported at the **first mapped base** of each alignment (the 5' end of
the alignment as seen by bamCoverage, using `--Offset 1 1`). Because reads carry a
1-base 5'-soft-clipped G, the actual G nucleotide sits **outside** the alignment:

- Forward strand: the G is 1 bp **upstream** of the bigWig position (position − 1)
- Reverse strand: the G is 1 bp **downstream** of the bigWig position (position + 1)

When extracting initiator dinucleotides in R (e.g. using rtracklayer + BSgenome),
shift the genomic coordinates by 1 in the appropriate direction before sequence lookup.

Example (R pseudocode):

```r
# Forward strand: TSS signal at position X means G is at X-1
fwd_tss <- resize(fwd_peaks, width = 1, fix = "start")
fwd_dinuc_range <- shift(fwd_tss, shift = -1)  # move to G position
fwd_dinuc_range <- resize(fwd_dinuc_range, width = 2)  # G and -1 base

# Reverse strand: TSS signal at position X means G is at X+1
rev_tss <- resize(rev_peaks, width = 1, fix = "end")
rev_dinuc_range <- shift(rev_tss, shift = 1)
rev_dinuc_range <- resize(rev_dinuc_range, width = 2)
```
