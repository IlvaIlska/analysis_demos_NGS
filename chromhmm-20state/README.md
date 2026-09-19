# ChromHMM 20-State Chromatin Annotation

Learns a 20-state chromatin annotation (ChromHMM) from a panel of 46 ChIP-seq marks in naïve
mouse B cells — histone PTMs, chromatin remodelers (Brg1, CHD4, MLL1, WDR5), architectural
proteins (CTCF, Rad21), HDAC1/2, p300, H2A.Z, and this portfolio's own 5hmC calls — then calls
per-bin chromatin states genome-wide. Source: `MH050_chromHMM_land/` (see repo root README for
provenance).

## Pipeline stages

| Stage | Script | What it does |
|---|---|---|
| 1 | `scripts/01_align_single_end.sh` | Trim (Trim Galore) → align to mm10 (Bowtie2, single-end) → filter to mapped reads → dedup (Picard) |
| 2 | `scripts/02_call_peaks.sh` | MACS2 `callpeak`, one mark's IP BAM vs. the shared input BAM, `-f BAM -q 0.05 --keep-dup all` |
| 3 | `scripts/03_binarize.sh` | ChromHMM `BinarizeBed -peaks`: turn all 46 marks' narrowPeaks into ChromHMM's binary presence/absence matrix, per chromosome |
| 4 | `scripts/04_learn_model.sh` | ChromHMM `LearnModel`: fit a 20-state HMM and produce the full report (segmentation, emission/transition parameters, TSS/TES enrichment, browser tracks, summary webpage) |

Each script has a header block listing its declared inputs, outputs, and parameters.

## No original driver script survived for stages 3–4 either

Same situation as the flagship pipeline's DiffBind stage: only the narrowPeak files, the
`cellmark.tab` sample table, and the `LearnModel` output (`model_20.txt`, `emissions_20.txt`,
etc.) survive from the original run — no `.sh` driver for the `BinarizeBed`/`LearnModel` calls
themselves. `scripts/03_binarize.sh` and `scripts/04_learn_model.sh` are written directly from
ChromHMM's own CLI usage text (`java -jar ChromHMM.jar BinarizeBed`/`LearnModel` with no
arguments prints full usage) plus the surviving `model_20.txt` naming (`nB`, 20 states, mm10),
not copied from any recovered script.

## Demo data

**Public accession:** [`GSE82144`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82144)
/ ENA study `SRP075985` (Kieffer-Kwon/Casellas mouse B-cell resource) — its `rB_wt_<mark>`
("resting B, wild-type" = naïve B) ChIP-seq runs cover **all 45 of the 46 marks** in this
pipeline's panel, plus one shared input control. Single-end, Illumina HiSeq 2500,
~25–30M reads/run. See `config/samples.tsv` for the full mark → run accession mapping (built from
that study's run-level metadata) and `config/cellmark_template.tab` for the ChromHMM-format
sample table. Where a mark has multiple replicates in the source (CTCF, Rad21, HDAC1), the demo
uses rep1; alternates are noted in `config/samples.tsv`.

The one mark public ChIP-seq data can't supply is `Son_5hmC` — this pipeline's only
self-generated mark. It comes directly from
[`../flagship-mnase-5hmc-diffbind/`](../flagship-mnase-5hmc-diffbind/) stage 03's output (the
sonicated-fraction 5hmC MACS2 peaks), so running this pipeline's demo end-to-end means running
the flagship's stages 0–3 first for at least one sonicated sample.

Fetch a mark's FASTQ with, e.g.:

```bash
prefetch SRR3619348 && fasterq-dump --split-files SRR3619348   # single-end: only *_1 is produced
```

## Environment

```bash
conda env create -f environment.yml
conda activate chromhmm-panel
```

ChromHMM itself is a jar, not a conda package (see the note in [`environment.yml`](environment.yml))
— point `$CHROMHMM_JAR` (or the `-j` flag on stages 3–4) at your own copy, e.g. download from
[compbio.mit.edu/ChromHMM](http://compbio.mit.edu/ChromHMM/). ChromHMM ships its own mm10 chrom
sizes/TSS/TES coordinate files, so `-a mm10` in stage 4 works without any extra genome download.
