# MNase-Fractionated 5hmC Profiling + DiffBind Differential Occupancy

CUT&Tag for 5-hydroxymethylcytosine (5hmC) performed on three chromatin fractions obtained by
partial micrococcal nuclease (MNase) digestion — **mono**nucleosome, **sub**mononucleosome
(sub-nucleosomal fragments), and **sonicated** (fragmented, MNase-independent) chromatin — from
naïve mouse B cells, three biological samples (`1590`, `1592`, `B6`). The pipeline aligns and
filters reads, calls peaks against matched input controls, and uses DiffBind to classify
5hmC-enriched regions by which chromatin fraction they're associated with.

This supports the analysis in Hsu et al., *BMC Genomics* 2026 (see repo root README for full
citation): the paper's central finding that most 5hmC sits on nucleosome-associated DNA.

## Pipeline stages

| Stage | Script | What it does |
|---|---|---|
| 0 | `scripts/00_fastqc.sh` | Raw-read FastQC |
| 1 | `scripts/01_trim_align_filter.sh` | Trim (Trim Galore) → align to mm10 (Bowtie2) → keep properly-paired reads → dedup (Picard) → RPKM bigWig |
| 2 | `scripts/02_t4_spikein_normalization.sh` | Align trimmed reads to a T4 phage spike-in genome → compute a scale factor from spike-in read depth → scaled bigWig for cross-sample normalization |
| 3 | `scripts/03_macs2_call_peaks.sh` | MACS2 `callpeak`, IP vs. matched input, `-f BAMPE -q 0.05 --keep-dup all` |
| 4 | `scripts/04_diffbind_differential_occupancy.R` | DiffBind: count reads in a consensus peakset, normalize, run all three pairwise fraction contrasts, classify sites as fraction-enriched (`posFold`/`negFold`) or non-significant (`NS`) at two fold thresholds |

Each script has a header block listing its declared inputs, outputs, and parameters — read that
before running rather than relying on this table alone.

## Stage 4 is a reconstruction, not a copy

The original DiffBind driver script for this analysis was never saved to disk — only a 504MB
`.RData` workspace and the downstream output directory structure survive
(`DiffBind_output/<date>/fold_1.5/` and `fold_2/`, each containing `posFold_*`, `negFold_*`, and
`NS_*` bed files per pairwise comparison). `scripts/04_diffbind_differential_occupancy.R` was
written from scratch to reproduce that structure, using:

- the sample sheet columns (`bamReads`/`bamControl`/`Peaks`, `Treatment` = mono/submono/sonicated,
  `Replicate` = 1590/1592/B6) — see `config/samples.tsv`
- the standard DiffBind workflow (`dba` → `dba.count` → `dba.normalize` → `dba.contrast` →
  `dba.analyze` → `dba.report`)
- the three pairwise comparisons implied by the surviving output (`monoVsubmono`, `monoVsonicated`,
  `sonicatedVsubmono`) and the two fold thresholds (1.5, 2) used to bin calls into `posFold`/`negFold`/`NS`

Treat this script's exact parameter defaults (FDR cutoff, normalization method) as reasonable
DiffBind defaults rather than a byte-for-byte recovery of what was originally run.

## Demo data

**Public accession:** [`SRP684244`](https://www.ncbi.nlm.nih.gov/sra/?term=SRP684244) — this
paper's own MNase-fraction 5hmC CUT&Tag data, *Mus musculus*, Illumina NovaSeq 6000.

12 runs are public, covering the **submono** and **sonicated** fractions (IP + input, ×3
biological samples each) — the **mono** fraction was not deposited. That means the public demo
can run the full stage 0–3 pipeline for all 12 runs, but stage 4 (DiffBind) demo is limited to
the `sonicatedVsubmono` contrast; `monoVsubmono` and `monoVsonicated` need the original mono-
fraction data to run for real (the script still runs against whatever fractions are present in
`config/samples.tsv`).

| Run accession | Library | Fraction | Sample | Type |
|---|---|---|---|---|
| SRR37650364 | submono_1590 | submono | 1590 | IP |
| SRR37650365 | submono_1590_input | submono | 1590 | input |
| SRR37650362 | submono_1592 | submono | 1592 | IP |
| SRR37650363 | submono_1592_input | submono | 1592 | input |
| SRR37650370 | submono_B6 | submono | B6 | IP |
| SRR37650371 | submono_B6_input | submono | B6 | input |
| SRR37650372 | sonicated_1590 | sonicated | 1590 | IP |
| SRR37650373 | sonicated_1590_input | sonicated | 1590 | input |
| SRR37650368 | sonicated_1592 | sonicated | 1592 | IP |
| SRR37650369 | sonicated_1592_input | sonicated | 1592 | input |
| SRR37650366 | sonicated_B6 | sonicated | B6 | IP |
| SRR37650367 | sonicated_B6_input | sonicated | B6 | input |

Fetch with, e.g.:

```bash
prefetch SRR37650364 && fasterq-dump --split-files SRR37650364
```

or pull directly over FTP from ENA (see `config/samples.tsv` for the file-report query used to
build this table: `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP684244&result=read_run&fields=run_accession,sample_title,library_name,fastq_ftp`).

References needed: a Bowtie2 index for `mm10`, and a Bowtie2 index for the T4 phage genome
(spike-in; any T4 phage reference FASTA works, e.g. RefSeq `NC_000866.4`).

## Environment

Original scripts assumed an OSC HPC environment (`module load fastqc bowtie2 samtools picard`,
a conda env for MACS2/deepTools, a separately-installed Trim Galore). For portability, the
cleaned scripts here document required tools by name in their header comments rather than
hardcoding `module load` — install via [`environment.yml`](environment.yml):

```bash
conda env create -f environment.yml
conda activate mnase-5hmc
```
