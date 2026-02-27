# TengLab-project

Reproduced code, pipelines, and documentation for two liquid biopsy papers studied in the Teng Lab.

---

## Subprojects

### [`songetal/`](songetal/README.md) — Gastric Carcinoma cfDNA Multi-Omic Pipeline

Code and data reference for:
> Song et al. "Plasma cfDNA multi-omic biomarkers profiling for detection and stratification of gastric carcinoma." *BMC Cancer* 25, 1003 (2025). https://doi.org/10.1186/s12885-025-14409-0

- **`songetal/scr/get_motif_fragsize.py`** — Extracts cfDNA fragmentation short/long ratios and 4-mer end-motif profiles from BAM files.
- **`songetal/scr/WGS.nf`** — Nextflow pipeline: QC (fastp) → alignment (BWA) → feature extraction.
- **`songetal/data/`** — Supplementary tables (WGS QC metadata, TF reference lists). Large feature matrix excluded from version control.

---

### [`ScanTecc/`](ScanTecc/README.md) — Cell-Free eccDNA Cancer Detection

Code and documentation for:
> Fang et al. "Detection of primary cancer types via fragment size selection in circulating cell-free extrachromosomal circular DNA." *Genome Medicine* 18, 18 (2026). https://doi.org/10.1186/s13073-025-01595-6

- **`ScanTecc/scantecc_eccDNA_identification_pipeline/`** — Multi-process eccDNA junction detection from BAM files (`efp_parallel.py`), CIGAR parsing (`parser.py`), and interval arithmetic (`genomic_interval.py`).
- **`ScanTecc/Figure*.ipynb` / `ScanTecc/Figure*.py`** — Notebooks and scripts reproducing all main figures (eccDNA abundance, size selection, gene annotation, ScanTecc classifier ROC curves).
- **`ScanTecc/Sfigure*.ipynb`** — Supplementary figure notebooks (RCA comparison, size distribution, gene density correlation).

---

## Data

Large data files and supplementary tables are **not tracked** in this repository (see `.gitignore`). See each subproject's `README.md` for download links and sources.
