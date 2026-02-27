# Song et al. (2025) — Gastric Carcinoma cfDNA Multi-Omic Pipeline

Documentation for the data and code reproduced from:

> Song et al. "Plasma cfDNA multi-omic biomarkers profiling for detection and stratification of gastric carcinoma." *BMC Cancer* 25, 1003 (2025). https://doi.org/10.1186/s12885-025-14409-0

---

## Project Structure

```
songetal/
├── data/
│   ├── Gastriccancerfeatures.csv   # Multi-omic feature matrix (733 samples × 5,644 features)
│   ├── MOESM2_ESM.csv              # Supplementary Table S1 — WGS mapping QC per sample
│   └── MOESM3_ESM.csv              # Supplementary Figure S1 — end-motif reference list
└── scr/
    ├── get_motif_fragsize.py       # Feature extraction script (motif + fragmentation)
    ├── WGS.nf                      # Nextflow pipeline for WGS pre-processing
    ├── nextflow.config             # Nextflow configuration
    └── explore.ipynb               # Exploratory analysis notebook
```

---

## Data Files

### `Gastriccancerfeatures.csv` — Multi-Omic Feature Matrix

The main machine-learning-ready feature table used to train the GC detection model.

| Property | Value |
|---|---|
| Rows (samples) | 733 |
| Columns (features) | 5,644 |
| Sample ID column | `ID` |
| Feature columns | `chr1_15000001_20000000` … `chr22_50000001_50500000` |

**Feature composition** (5,643 feature columns):
- **Fragmentation profile** (472 features): Short/long fragment ratios in 5 Mb genomic bins across the mappable human genome (hg38).
- **End-motif profile** (256 features): Proportional abundance of all 256 possible 4-mer sequences at cfDNA fragment R1 ends. In the paper, 44 representative motif codes are highlighted (see `MOESM3_ESM.csv`).
- **Copy Number Variation (CNV)** (4,915 features): Log2 CNV ratios in non-overlapping 0.5 Mb bins genome-wide, computed with ichorCNA.

Together these 5,643 features (472 + 256 + 4,915) are used in the XGBoost classifier described in the paper, with feature weights of ~0.305 (fragmentation), ~0.329 (end-motif), and ~0.366 (CNV).

---

### `MOESM2_ESM.csv` — Supplementary Table S1: WGS Mapping QC

Per-sample sequencing quality control and cohort metadata for all 733 participants.

| Column | Description |
|---|---|
| `Samples` | Sample identifier |
| `Raw reads` | Total paired-end reads before QC |
| `Raw bases` | Total base pairs sequenced |
| `clean reads` | Reads passing QC (adapter trimmed, quality filtered) |
| `Q20(%)` | Bases with Phred score ≥ 20 (count + %) |
| `Q30(%)` | Bases with Phred score ≥ 30 (count + %) |
| `Clean bases` | Total clean base pairs |
| `clean rate(%)` | Fraction of reads retained after QC |
| `map-rate` | Fraction of reads mapped to hg38 |
| `depth` | Mean genome coverage (average ~1.1×) |
| `cohort` | `healthy`, `GC`, or `Bengin` |
| `pathology` | Specific diagnosis for non-healthy samples |
| `Stage` | Pathological stage for GC samples (I–IV) |

**Cohort breakdown:**

| Cohort | n |
|---|---|
| Healthy | 521 |
| Gastric Cancer (GC) | 131 |
| Benign gastric disease | 81 |
| **Total** | **733** |

**GC staging (52 of 131 patients with staging information):**

| Stage | n |
|---|---|
| I | 1 |
| II | 15 |
| III | 16 |
| IV | 20 |

**Benign disease types** include: gastric polyps, gastrointestinal polyps, gastric ulcer, gastroduodenal ulcers, esophagitis, colorectal polyps, GIST.

---

### `MOESM3_ESM.csv` — Supplementary Figure S1: End-Motif Reference

A two-column reference table listing transcription factors relevant to the GC analysis.

| Column | Description |
|---|---|
| `AnimalTFDB4 (Homo sapiens)` | Full set of human TFs from the AnimalTFDB4 database that were differentially expressed between healthy and GC TCGA-STAD samples (128 entries) |
| `GC-specific TFs` | Subset of 15 TFs identified as GC-specific based on TFBS coverage differences across cohorts |

**GC-specific TFs** (column 2): `ARID3A`, `CTCFL`, `EGR1`, `EGR2`, `ETV4`, `FEZF1`, `FOXM1`, `HLF`, `KLF15`, `SALL4`, `TBX5`, `WT1`, `ZBTB16`, `ZIC2`

Note: `CTCFL`, `EGR1`, `EGR2`, `ETV4`, `FOXM1`, and `KLF15` were previously reported as GC-associated; `FEZF1`, `TBX5`, and `WT1` are novel findings from this study.

---

## Code

### `get_motif_fragsize.py` — cfDNA Feature Extraction

**Purpose:** Extracts two cfDNA biomarker features from a BAM alignment file — the fragmentation short/long ratio profile and the 4-mer end-motif profile — for each genomic interval provided.

**Usage:**
```bash
python get_motif_fragsize.py <input.bam> <output_prefix> <genomic_intervals.csv>
```

| Argument | Description |
|---|---|
| `sys.argv[1]` | Path to sorted, indexed BAM file |
| `sys.argv[2]` | Output file prefix |
| `sys.argv[3]` | CSV file of genomic intervals (columns: `chr`, `start`, `end`) |

**Read filtering:**
- Flag mask `3332` excludes unmapped, secondary, supplementary, QC-failed, and duplicate reads.
- Mapping quality > 30 required.
- Only read 1 (`is_read1`) with a non-zero insert size (`isize`) is processed.

**Fragment size classification:**
- Short fragments: insert size 100–150 bp (nucleosome-free, enriched in ctDNA)
- Long fragments: insert size 151–220 bp (nucleosome-wrapped)

**End-motif extraction:**
- For read 1: first 4 bases of `query_sequence` (5′ end motif).
- 4-mers containing `N` are discarded.
- All 4-mer counts are normalized to proportions (sum = 1).

**Outputs** (four TSV files):

| File | Content |
|---|---|
| `<prefix>_motif_original.csv` | Raw motif proportions — columns: `tag` (4-mer), `num` (proportion) |
| `<prefix>_motif_normalize.csv` | Mean-centered motif proportions (subtract genome-wide mean per motif) |
| `<prefix>_fragment_original.csv` | Raw short/long ratios per interval — columns: `chr`, `start`, `end`, `Ratio` |
| `<prefix>_fragment_normalize.csv` | Mean-centered short/long ratios (subtract per-sample mean across bins) |

**Normalization method:** Both motif and fragment features are normalized by subtracting the per-sample mean across all intervals/motifs (mean-centering). This corrects for global depth differences between samples.

---

### `WGS.nf` — Nextflow Pre-Processing Pipeline

**Purpose:** End-to-end Nextflow (DSL1) workflow that takes raw FASTQ files and produces sorted, duplicate-marked BAM files, then runs motif/fragmentation feature extraction.

**Pipeline stages:**

1. **QC** (`fastp`): Adapter trimming and quality filtering of paired-end reads. Applies a minimum length filter and outputs clean FASTQ files and QC statistics.
2. **check_QC**: Validates that the clean read rate passes a configurable threshold.
3. **Map** (`BWA mem` + `sambamba`): Aligns clean reads to hg38, sorts, marks PCR duplicates, merges per-chromosome BAMs, and indexes the final BAM. Computes mapping rate with `samtools flagstat`.
4. **Motif**: Calls `get_motif_fragsize.py` on the final BAM to produce `*motif_original.csv` and `*fragment_original.csv` per sample.

**Key tools used:** `fastp`, `BWA mem`, `sambamba`, `samtools`, `sinotools` (duplicate marking)

**Input:** Paired FASTQ files matching `*{_1,_2}*{fq,fastq}*` in the specified input directory.

**Output per sample:** Sorted/marked BAM + index, mapping QC stats, motif CSV, fragmentation ratio CSV.

---

## Biological Context

This repository reproduces features from a non-invasive liquid biopsy study for gastric carcinoma (GC) detection. Plasma cell-free DNA (cfDNA) from 733 participants was sequenced at low depth (~1.1×) and three complementary feature types were extracted:

1. **Fragmentation profile** — cfDNA from tumors tends to be shorter than cfDNA from normal cells, reflecting differential nucleosome positioning and increased nuclease activity. Short/long ratios across 5 Mb bins reveal genome-wide fragmentation patterns distinct in GC.

2. **End-motif profile** — The 4-mer sequence at each cfDNA fragment end reflects the sequence preferences of the nucleases that generated it, which differ between healthy and tumor tissue environments.

3. **CNV profile** — Tumor cells accumulate chromosomal copy number gains and losses. Even at 1.1× depth, these are detectable genome-wide and are highly discriminatory for GC (relative feature weight ~0.366).

An XGBoost classifier trained on these 5,643 features achieves **AUC = 0.998**, sensitivity **94.87%**, and specificity **99.35%** for distinguishing healthy individuals from GC patients.
