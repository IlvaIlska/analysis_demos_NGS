# ChromHMM 20-State Chromatin Annotation

**Status: scaffolded, not yet built out.** This is one of two secondary pipelines held in reserve
(see repo root README) — organized here so the directory structure and source material are
locked in, but the stage scripts themselves haven't been cleaned up yet the way the
[flagship pipeline](../flagship-mnase-5hmc-diffbind/) has.

## What it does

Learns a multi-state chromatin annotation (ChromHMM) from a panel of ~43 ChIP-seq/CUT&Tag marks
in naïve mouse B cells — histone PTMs, chromatin remodelers (Brg1, CHD4, MLL1, WDR5), architectural
proteins (CTCF, Rad21), HDAC1/2, p300, and H2A.Z — then calls per-mark, per-state genome
annotations from the learned model. Source: `MH050_chromHMM_land/` (see below).

## Pipeline stages (from source material)

1. **Binarize** each mark's narrowPeak calls against the genome, keyed by a `cellmark.tab`
   sample table (`<celltype> <mark> <narrowPeak file>`, one row per mark — e.g.
   `062323_cellmark.tab` lists 43 marks for cell type `nB`). ChromHMM's `BinarizeBed`.
2. **LearnModel**: fit a hidden Markov model over the binarized marks to learn emission/transition
   parameters and assign chromatin states (20 states in the original run — `062323_learningOutput_all/`).
3. State calling / annotation: apply the learned model to produce a per-bin state assignment
   across the genome, and characterize what each state's mark combination represents (active
   promoter, enhancer, heterochromatin, etc.).

Variant runs also exist that drop specific marks to test robustness (`061223_binarized_MH049_hmCless`,
`062323_binarized_MH049_H2AZless`) — worth keeping as a "does the model degrade gracefully"
talking point if this pipeline gets built out further.

## Demo data (candidate public accessions — not yet wired up)

- [`GSE116208`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116208) — "TET enzymes
  augment AID expression via 5hmC modifications at the Aicda superenhancer" (Lio lab). Likely
  source for some of the ChromHMM mark panel.
- [`GSE82144`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82144) — "Myc regulates
  chromatin decompaction and nuclear architecture during B cell activation"
  (Kieffer-Kwon/Casellas). Source of CTCF, Rad21, Brg1, MLL1, WDR5, HDAC1/2, GCN5, CHD4, p300
  peaks present in `MH050_chromHMM_land/` as `*_0h.narrowPeak`.

A full 43-mark demo isn't realistic to assemble from public data alone; a reduced-panel demo
(the subset actually traceable to GSE116208/GSE82144) is the more honest target once this
pipeline gets built out.

## TODO before this is portfolio-ready

- [ ] Identify which of the 43 marks in `062323_cellmark.tab` map to which public accession/run
- [ ] Write a cleaned `01_binarize.sh` (ChromHMM `BinarizeBed`, parameterized on a cellmark table + peaks dir)
- [ ] Write a cleaned `02_learn_model.sh` (ChromHMM `LearnModel`, parameterized on state count)
- [ ] Write a `03_annotate_states.sh` / state-characterization step
- [ ] Decide whether to keep the full 20-state model or trim to states with clear public-data support
