PROJECT: cen-tss-analysis
AUTHOR:  Pavel Nikitin (pnikitin, based on old script paths)
PURPOSE: Identify and characterise transcription start sites (TSS) in human
         centromeric alpha-satellite (HOR) regions using CAGE data remapped
         to the T2T CHM13 reference genome. Part of a thesis project.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
SCIENTIFIC BACKGROUND
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

- CAGE (Cap Analysis of Gene Expression) captures the 5' end of capped RNA,
  identifying TSSs at single-nucleotide resolution.
- CAGE reads often start with a 5'-soft-clipped G because the CAGE template-
  switching mechanism adds a G to the 5' end of the cDNA. This G is the
  initiator nucleotide (+1 position of the TSS).
- The pipeline remaps CAGE reads to a specialised centromeric reference
  (t2t_hor_with_strand_centromeres_ref) using bowtie2 local mode, retaining
  up to 800 multimapping alignments (-k 800) to handle repetitive sequences.
- Cell lines/samples processed so far (from scripts_old/): H9 hESC, K562,
  HepG2, HeLaS3, GM12878, MCF7, LUHMES (neural differentiation time series,
  day 1/3/6, CAGE and NETCAGE), and others. All are human cell lines.
- The centromeric reference contains HOR (higher-order repeat) units with
  strand information. The transposon BED file marks centromeric TEs to exclude.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PIPELINE OVERVIEW (6 steps, bash, no workflow manager)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Input:  STAR-mapped CAGE BAM (full genome, T2T CHM13).
        The pipeline begins here. Upstream steps (trimming with Trim Galore,
        STAR mapping via CAGEseq2/Nextflow) are handled separately.

Step 1  get_remap.sh          Extract reads in centromeric BED regions + all
                              unmapped reads → merged FASTQ.gz for remapping.
Step 2  bowtie2_remap.sh      Remap to centromeric reference, bowtie2 local
                              mode, -k 800, multimappers retained.
Step 3  select_tss_alignments.sh  Select reads with exactly 1 5'-soft-clipped
                              G (pyselectal.py -n 1 -m 1 -x G). Checks both
                              CIGAR pattern AND actual nucleotide identity.
Step 4  xm_filter.sh          Remove reads with >= 10% mismatches (XM tag /
                              read length >= 0.10).
Step 5  filter_transposons.sh Remove reads overlapping centromeric transposons
                              (bedtools intersect -v -abam).
Step 6  bamcoverage_5p.sh     Produce SAMPLE.fwd.bw and SAMPLE.rev.bw:
                              strand-specific 5'-end bigWigs, 1 bp bins,
                              --Offset 1 1, both files use positive values.

Output: SAMPLE.fwd.bw, SAMPLE.rev.bw, SAMPLE.log (per sample).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
REPOSITORY LAYOUT
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

run_pipeline.sh               Main orchestration script (see interface below)
pipeline.conf                 All configurable paths and parameters
scripts/get_remap.sh          Step 1
scripts/bowtie2_remap.sh      Step 2
scripts/select_tss_alignments.sh  Step 3 (wrapper for pyselectal.py)
scripts/pyselectal.py         Python BAM filter tool (checks CIGAR + base)
scripts/xm_filter.sh          Step 4
scripts/filter_transposons.sh Step 5
scripts/bamcoverage_5p.sh     Step 6 (produces .fwd.bw and .rev.bw)
scripts_old/                  Original scripts, kept for reference only
README.md                     Full documentation

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
RUN INTERFACE
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

bash run_pipeline.sh -i INPUT -o OUTPUT_DIR [options]

  -i / --input PATH          BAM file or directory of BAM files
  -o / --output DIR          Output root directory
  -c / --config FILE         Config file (default: pipeline.conf)
  -t / --threads N           Thread count (overrides THREADS in config)
  --keep-intermediates       Do not delete intermediate BAM/FASTQ files
  --resume                   Skip steps whose output already exists
  --dry-run                  Print commands, do not execute

Examples:
  bash run_pipeline.sh -i data/bams/ -o results/ -c pipeline.conf -t 16
  bash run_pipeline.sh -i data/sampleA.bam -o results/ --resume
  bash run_pipeline.sh -i data/bams/ -o results/ --dry-run

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
PIPELINE.CONF — MUST FILL IN BEFORE RUNNING
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

CENTROMERE_BED        censat_hor_with_strand_R.bed           (absolute path)
TRANSPOSON_BED        adjusted_centromeric_transposons_full_ranges.bed (absolute path)
BOWTIE2_INDEX_PREFIX  t2t_hor_with_strand_centromeres_ref    (prefix, no extension)
THREADS               20 (default)

BOWTIE2_EXTRA_ARGS and PYSELECTAL_ARGS are bash arrays — edit as arrays,
not as quoted strings, otherwise word-splitting breaks the argument passing.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
DOWNSTREAM ANALYSIS PLANS (R)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

After running the pipeline, the bigWig files will be loaded in R for:
- TSS peak calling / clustering (likely CAGEr or similar)
- Initiator dinucleotide extraction: the nucleotide at the TSS (+1, the G)
  and the nucleotide immediately upstream (-1 position).

IMPORTANT OFFSET NOTE:
  The bigWig coverage is at the FIRST MAPPED BASE (after the soft-clip),
  not at the G itself. For reads with 1S soft-clip:
    forward strand: G is 1 bp UPSTREAM  of the bigWig position (X - 1)
    reverse strand: G is 1 bp DOWNSTREAM of the bigWig position (X + 1)
  Shift coordinates accordingly in R before sequence extraction (rtracklayer
  + BSgenome). See README section 5 for R pseudocode.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
KEY DESIGN DECISIONS (do not change without good reason)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. Both centromeric reads AND unmapped reads go into the bowtie2 remapping.
   This is intentional — needed to capture reads that STAR missed.

2. 5'-G filter checks CIGAR + actual nucleotide (pyselectal.py -n 1 -m 1 -x G).
   The original scripts only checked CIGAR. This is a deliberate improvement.

3. Intermediate files are deleted after each downstream step succeeds, NOT
   via trap. Failed runs leave the last successfully created file on disk.

4. process_sample() runs in a subshell; one failing sample does not abort
   the batch. All failures reported at the end with exit code 1.

5. filter_transposons uses bedtools intersect -v (not grep by read name).

6. BigWigs are two separate files with positive values (not combined/negative).

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CURRENT STATUS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

- All scripts written and refactored. NOT yet tested on real data.
- pipeline.conf paths are placeholders — must be updated before first run.
- pysam Python package must be installed (pip install pysam or conda).
- The bowtie2 index must exist at BOWTIE2_INDEX_PREFIX.
- Next step: set paths in pipeline.conf and do a --dry-run, then a real run.

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
ORIGINAL SCRIPTS SUMMARY (scripts_old/ — reference only, not executed)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

bowtie2_remap.sh        Hardcoded sample list + paths. Replaced by step 2.
get5GsoftclipBAM.sh     CIGAR-only soft-clip filter (no G base check). Step 3 improves this.
get5GsoftclipBAM_initial.sh  Same filter applied to STAR BAM. Not in this pipeline.
xm_filter.sh            Active XM filter (< 10% mismatches). Much of the body is
                        commented-out code with a stricter G/C base check.
star.sh                 STAR mapping — done upstream with Nextflow, not here.
trimgalore.sh           Adapter trimming — done upstream with Nextflow, not here.
wget.sh                 Data download scripts for public datasets (SRA/ENCODE).
get5Greads.sh           FASTQ-level pre-mapping G filter. Older/abandoned approach.
