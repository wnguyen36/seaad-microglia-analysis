# SEA-AD Microglia scRNA-seq Analysis

Single-cell transcriptomic analysis of **240,651 microglia** from the Seattle Alzheimer's Disease Brain Cell Atlas (SEA-AD), investigating disease-associated microglial (DAM) activation states across Alzheimer's disease severity, APOE genotype, and neuropathological burden.

---

## Overview

Microglia are the brain's resident immune cells and a central player in Alzheimer's disease (AD) pathology. This project applies a complete single-cell RNA-seq analysis pipeline to the SEA-AD microglia dataset — one of the largest, most deeply annotated human brain single-cell datasets available — to characterize how microglial cell states shift with disease progression and genetic risk.

**Key questions addressed:**
- Do disease-associated microglia (DAM) accumulate with increasing Braak stage and cognitive decline?
- Do APOE ε4 carriers show elevated DAM activation at the donor level?
- Which microglial supertypes drive the transition from homeostatic to disease-associated states?
- What genes are differentially expressed between DAM-high and DAM-low donors?

---

## Dataset

| Property | Value |
|---|---|
| Source | Allen Institute for Brain Science — SEA-AD consortium |
| Cell type | Microglia / perivascular macrophages |
| Cells | 240,651 |
| Genes | 36,601 |
| Donors | 84 |
| Normalization | Log-normalized (log1p), pre-processed by Allen Institute |

**Clinical metadata integrated:**
- Cognitive status (No dementia / Dementia)
- Braak neurofibrillary tangle stage (0 – VI)
- CERAD neuritic plaque score (Absent → Frequent)
- Thal amyloid phase (0 – 5)
- Overall AD neuropathological change (Not AD → High)
- APOE genotype (2/2, 2/3, 3/3, 2/4, 3/4, 4/4)
- Severely affected donor flag
- Sex, age at death, post-mortem interval (PMI)
- Quantitative neuropathology: Iba1+ area, AT8+ tau burden, 6E10+ amyloid, NeuN+ cell density, RIPA pTau/Aβ42

**Microglial supertypes (Allen Institute annotations):**

| Supertype | Cells | Description |
|---|---|---|
| Micro-PVM_2 | 141,248 | Canonical homeostatic microglia |
| Micro-PVM_2_3-SEAAD | 53,815 | SEAAD-enriched transitional subtype |
| Micro-PVM_3-SEAAD | 29,805 | SEAAD-enriched activated subtype |
| Micro-PVM_1 | 7,909 | Canonical homeostatic microglia |
| Lymphocyte | 4,607 | Contaminating lymphocytes |
| Micro-PVM_4-SEAAD | 1,767 | SEAAD-enriched reactive subtype |
| Micro-PVM_2_1-SEAAD | 1,458 | SEAAD-enriched subtype |
| Monocyte | 42 | Contaminating monocytes |

---

## Pipeline

The analysis is structured as a modular Python package (`src/`) with a primary analysis notebook (`notebooks/SeaAD_Analysis.ipynb`).

```
SEA-AD-scRNAseq/
├── src/
│   ├── preprocess.py     # QC filtering, normalization, checkpointing
│   ├── cluster.py        # PCA, KNN graph, UMAP, Leiden clustering
│   ├── annotate.py       # Marker gene identification, cell type assignment, gene set scoring
│   └── visualize.py      # QC, HVG, PCA, UMAP, dotplot visualizations
├── notebooks/
│   ├── SeaAD_Analysis.ipynb   # Main analysis (Cells 1–12)
│   ├── 01_qualcontrol.ipynb
│   ├── 02_normalization.ipynb
│   ├── 03_clustering.ipynb
│   └── 04_labeling.ipynb
├── data/
│   ├── SEAAD_microglia.h5ad        # Primary dataset (Allen Institute)
│   ├── SEAAD_donor_metadata.xlsx   # Clinical donor metadata
│   └── SEAAD_neuropathology.csv    # Quantitative neuropathology measurements
├── figures/                        # Publication-quality figures (300 dpi)
└── results/                        # Saved AnnData objects and CSVs
```

### Analysis steps

| Cell | Step | Key tools |
|---|---|---|
| 1 | Environment setup & module imports | Scanpy, pandas, NumPy, SciPy, seaborn |
| 2 | Dataset loading & metadata inspection | `sc.read_h5ad`, pandas |
| 3 | Clinical variable audit | Value counts across all AD staging variables |
| 4 | Allen Institute cell type annotation review | Subclass / Supertype / Class labels |
| 5 | Normalization status verification | Sparse matrix diagnostics, known marker validation |
| 6 | UMAP visualization by clinical variables | Supertype, cognitive status, Braak, severely affected donors |
| 7 | Neuropathology metadata integration | Donor-level join of Iba1, AT8, 6E10, NeuN, pTau/Aβ42 |
| 8 | Gene set scoring (DAM & Homeostatic) | `sc.tl.score_genes`, 27-gene DAM signature, 8-gene homeostatic signature |
| 8B | Supertype activation profiling | UMAP overlays, donor-aggregated bar charts, scatter plot |
| 8C | Statistical testing (supertype-level) | Mann-Whitney U, Benjamini-Hochberg FDR |
| 9 | Subtype composition analysis by Braak | Donor-level proportions, seaborn boxplots |
| 10 | Donor-level correlation analysis | Spearman ρ vs Braak, CERAD, Thal, ADNC, cognitive status, age |
| 11 | APOE4 subgroup analysis | Carrier vs non-carrier Mann-Whitney U at donor level |
| 12 | Pseudobulk differential expression | CPM normalization, Mann-Whitney per gene, BH correction, volcano plot |

---

## Gene Signatures

### Disease-Associated Microglia (DAM) — 27 genes
Curated from Keren-Shaul et al. (2017, *Cell*), Krasemann et al. (2017), and Haage et al. (2024). Mouse-to-human ortholog conversion performed via MGI records; mouse-specific genes without confirmed 1:1 human orthologs excluded.

**DAM stage 1:** `APOE`, `B2M`, `FTH1`, `CSTB`, `LYZ`, `CTSB`, `TYROBP`, `TIMP2`, `CTSD`

**DAM stage 2:** `CD9`, `CD63`, `SERPINE2`, `SPP1`, `CADM1`, `CD68`, `CTSZ`, `AXL`, `CLEC7A`, `CTSA`, `CD52`, `CSF1`, `LPL`, `CTSL`, `CST7`, `ITGAX`, `GUSB`, `HIF1A`

All 27 genes confirmed present in the 36,601-gene dataset.

### Homeostatic Microglia — 8 genes
`SALL1`, `HEXB`, `CX3CR1`, `TMEM119`, `TREM2`, `P2RY12`, `MERTK`, `PROS1`

Mouse-specific genes (`SIGLECH`, `GPR43`/`FFAR2`) excluded where no validated human ortholog exists.

---

## Key Results

### DAM Score Correlates with AD Neuropathological Severity

Donor-level Spearman correlations between mean DAM score and neuropathological variables (n = 84 donors):

| Variable | Spearman ρ | Significant |
|---|---|---|
| Braak Stage | Positive | Yes |
| CERAD Score | Positive | Yes |
| Thal Phase | Positive | Yes |
| Overall AD Neuropathological Change | Positive | Yes |
| Cognitive Status | Positive | Yes |
| Age at Death | — | — |

DAM scores increase monotonically across Braak stages 0 → VI, consistent with progressive microglial activation tracking tau pathology burden.

### SEAAD Subtypes Are More Activated Than Canonical Subtypes

Donor-aggregated DAM scores compared across supertypes using Mann-Whitney U with BH correction. SEAAD-enriched subtypes (`Micro-PVM_2_3-SEAAD`, `Micro-PVM_3-SEAAD`, `Micro-PVM_4-SEAAD`) show significantly elevated DAM scores and reciprocally reduced homeostatic scores relative to canonical `Micro-PVM_1` and `Micro-PVM_2` populations.

The supertype activation scatter plot shows a clear inverse DAM–Homeostatic axis, with SEAAD subtypes occupying the high-DAM / low-homeostatic quadrant.

### APOE4 Carriers Show Elevated DAM Activation

At the donor level (treating each of 84 donors as one independent observation):
- **APOE4 carriers** (3/4, 4/4, 2/4 genotypes)
- **APOE4 non-carriers** (3/3, 2/3, 2/2 genotypes)

Mann-Whitney U test on donor-level mean DAM scores demonstrates elevated activation in APOE4 carriers, consistent with the known role of APOE4 as the strongest genetic risk factor for late-onset AD.

### Pseudobulk Differential Expression: DAM-high vs DAM-low Donors

Pseudobulk analysis using raw UMI counts summed per donor, CPM-normalized, log1p-transformed, and tested per gene with Mann-Whitney U (BH-corrected):

| Direction | Genes | Threshold |
|---|---|---|
| Upregulated (DAM-high) | 27 | FDR < 0.05, mean diff > 0.5 log1p-CPM |
| Downregulated (DAM-high) | 22 | FDR < 0.05, mean diff < −0.5 log1p-CPM |

**Top upregulated genes:** `S100A6`, `LINC00482`, `TIMP1`, `LGALS1` — complement and inflammation-associated

Upregulated DAM signature genes overlap with known AD-risk and lysosomal activation programs. Results are exploratory (DAM grouping derived from the same expression data); held-out validation would be required for confirmatory claims.

---

## Figures

| Figure | Description |
|---|---|
| `umap_supertype.png` | UMAP colored by Allen Institute microglial supertypes |
| `umap_cognitive_status.png` | UMAP colored by dementia vs no dementia |
| `umap_braak_stage.png` | UMAP colored by continuous Braak stage |
| `umap_supertype_vs_scores.png` | Side-by-side UMAP: supertypes / DAM score / homeostatic score |
| `supertype_activation_profiles.png` | Donor-aggregated DAM + homeostatic bar charts and scatter |
| `composition_braak.png` | SEAAD subtype abundance by Braak stage (boxplots) |
| `DAM_correlations_donor_level.png` | Spearman ρ bar chart: DAM score vs neuropathological variables |
| `DAM_by_Braak_donor_level.png` | DAM score by Braak stage (donor-level boxplot) |
| `DAM_by_APOE_genotype_donor_level.png` | DAM score across all 6 APOE genotypes with per-group n |
| `APOE4_scores_donor_level.png` | DAM and homeostatic scores: APOE4 carriers vs non-carriers |
| `volcano_DE.png` | Volcano plot: pseudobulk DE, DAM-high vs DAM-low donors |

---

## Methods Summary

**Normalization:** Data provided by the Allen Institute is already log1p-normalized from raw UMIs. Confirmed via matrix value range (0–5.07) and `log1p` key in `adata.uns`.

**Gene set scoring:** `sc.tl.score_genes` (Seurat-style control-gene subtraction). Scores are computed per cell then aggregated to donor level by mean for all downstream statistical tests.

**Donor-level aggregation:** All statistical comparisons (supertype testing, Braak correlation, APOE analysis) operate on per-donor mean scores rather than pooled cells, preventing pseudoreplication inflating sample size.

**Pseudobulk DE:** Raw UMI counts (`adata.layers['UMIs']`) summed per donor → CPM → log1p. Mann-Whitney U per gene, Benjamini-Hochberg FDR correction. Thresholds: FDR < 0.05 and |mean difference| > 0.5.

**Multiple testing correction:** Benjamini-Hochberg FDR applied for all multi-gene and multi-group comparisons (`statsmodels.stats.multitest.multipletests`).

---

## Technologies

| Category | Tools |
|---|---|
| Core analysis | [Scanpy](https://scanpy.readthedocs.io/) 1.x, AnnData |
| Numerical | NumPy, SciPy (Spearman, Mann-Whitney U) |
| Data manipulation | pandas |
| Visualization | Matplotlib, seaborn |
| Statistics | SciPy stats, statsmodels (BH correction) |
| Data format | HDF5 / h5ad |
| Language | Python 3.10 |

---

## References

1. Keren-Shaul H, et al. (2017). A unique microglia type associated with restricting development of Alzheimer's disease. *Cell*, 169(7), 1276–1290. https://doi.org/10.1016/j.cell.2017.05.018
2. Krasemann S, et al. (2017). The TREM2-APOE pathway drives the transcriptional phenotype of dysfunctional microglia in neurodegenerative diseases. *Immunity*, 47(3), 566–581.
3. Haage V, et al. (2024). DAM revisited: new insights into microglial states in neurodegeneration.
4. Butovsky O, et al. (2014). Identification of a unique TGF-β–dependent molecular and functional signature in microglia. *Nature Neuroscience*, 17, 131–143.
5. Bennett ML, et al. (2016). New tools for studying microglia in the mouse and human CNS. *PNAS*, 113(12), E1738–E1746.
6. Lambert JC, et al. (2013). Meta-analysis of 74,046 individuals identifies 11 new susceptibility loci for Alzheimer's disease. *Nature Genetics*, 45, 1452–1458.
7. SEA-AD Consortium. Seattle Alzheimer's Disease Brain Cell Atlas. Allen Institute for Brain Science. https://portal.brain-map.org/explore/seattle-alzheimers-disease

---

## Data Access

The SEA-AD dataset is publicly available through the Allen Brain Cell Atlas portal. Donor metadata and neuropathology tables are distributed alongside the atlas.

> **Note:** The `.h5ad` file (`SEAAD_microglia.h5ad`, ~2 GB) is excluded from version control via `.gitignore`. Download from the Allen Brain Cell Atlas portal to reproduce this analysis.
