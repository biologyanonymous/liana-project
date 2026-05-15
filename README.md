# Tree Architecture × Liana Infestation — Analysis Code and Data

**Project:** Tree Architectural Feedbacks in Response to Liana Canopy Infestation  
**Study site:** 50-ha Forest Dynamics Plot, Barro Colorado Island, Panama  
**Code DOI:** *(Zenodo DOI will appear here after deposit)*

---

## Repository structure

```
mbiol-project/
├── data/                        ← analysis-ready CSV files (see below)
├── models/                      ← pre-fitted brms .rds files
├── figures/                     ← output figures (generated on render)
├── tables/                      ← output tables (generated on render)
├── analysis_and_figures.qmd     ← main analysis and figures
├── data_processing.qmd          ← data cleaning and processing pipeline
├── mm_waic_comparison.R         ← Michaelis–Menten vs log-log WAIC comparison
└── README.md                    ← this file
```

---

## Reproducing the analysis

Open `analysis_and_figures.qmd` in RStudio and click **Render**. All figures are written to `figures/` and tables to `tables/` automatically. Pre-fitted Bayesian models are stored as `.rds` files in `models/` and loaded directly; delete any `.rds` file and re-render to refit from scratch.

`data_processing.qmd` documents the upstream processing steps (field data cleaning, infestation-history classification, TLS-to-stem-map matching, crown metric extraction, PAI estimation). It is set to `eval: false` as the upstream inputs (raw TLS point clouds, dendrometer `.rdata` timeseries) are large files not included in this repository. Raw data are available on reasonable request.

### R version and key packages

| Package | Version | Purpose |
|---------|---------|---------|
| R | 4.5.2 | — |
| brms | ≥ 2.21 | Bayesian multilevel models |
| tidyverse | ≥ 2.0 | Data wrangling and plotting |
| ggplot2 | ≥ 3.5 | Figures |
| knitr / kableExtra | — | Tables |
| writexl | — | Excel output |

---

## Data files

### `dat_h1.csv` — H1 structural legacies
**n = 251 trees | 67 species | 24 subplots**

TLS crown metrics (2019) merged with infestation histories reconstructed from annual dendrometer records (2011–2019). Trees assigned to one of three mutually exclusive infestation-history categories.

| Column | Description |
|--------|-------------|
| `tag6` | 6-digit ForestGEO stem tag |
| `species` | Species code |
| `sbpltnm` | 40 × 40 m subplot ID |
| `H1_category` | `Never infested` / `Lost lianas` / `Persistently infested` |
| `dbh_cm` | Diameter at breast height (cm) |
| `height` | Tree height (m; TLS 2019) |
| `crown_depth` | Crown depth (m) |
| `pa` | Projected crown area (m²) |
| `crown_vol` | Crown volume (m³) |
| `slenderness` | Height / DBH |
| `cbh_ratio` | Crown base height / total height |
| `log_dbh`, `log_height`, `log_depth`, `log_area`, `log_vol` | Log-transformed metrics |

---

### `clean_pai.csv` — H1 plant area index
**n = 575 trees | 91 species | 24 subplots**

Voxel-based PAI estimated separately for host-tree points and combined host + liana points.

| Column | Description |
|--------|-------------|
| `tag` | Tree tag |
| `Species` | Species code |
| `sbpltnm` | Subplot ID |
| `dbh_cm` | DBH (cm) |
| `infest_cat` | `None` / `Light` (1–50%) / `Heavy` (51–100%) |
| `PAI_tree_only` | Host-only PAI (m² m⁻²) |
| `PAI_combined` | Combined host + liana PAI (m² m⁻²) |

---

### `dat_h2.csv` — H2 colonisation predictors (TLS)
**n = 181 trees | 61 species | 24 subplots**

Trees uninfested in 2019, classified by whether they gained liana infestation by 2025. Crown metrics are TLS 2019 measurements (pre-outcome baseline). Column structure identical to `dat_h1.csv` plus `H2_category` (`remained uninfested` / `gained infestation`).

---

### `dat_h2_2011.csv` — H2 colonisation predictors (2011 heights)
**n = 414 trees | 85 species**

Trees uninfested in 2011, classified by whether they gained liana infestation by 2019. Ground-measured heights from 2011 laser rangefinder surveys.

---

### `dat_mode.csv` — climbing mode × crown architecture
Fieldwork data (2025) joined with TLS crown metrics. One row per tree. Includes dominant climbing mode and active/passive climbing class.

---

### `pathway_dat.csv` — infestation pathway × colonisation outcome
Infestation pathway type (Direct / Lateral / Both) joined with H2 outcome. Used for chi-square test.

---

### `modes_clean.csv` — climbing mode records (long format)
One row per tree × observed liana climbing mode. Used for species-level proportional bar charts.

---

## Data provenance

| Source | Description |
|--------|-------------|
| TLS 2019 | Point clouds from 24 BCI 50-ha plot subplots; processed with RayCloudTools (Lowe & Stepanas, 2021) and ITSMe (Terryn et al., 2023) |
| Dendrometer records 2011–2025 | Annual liana infestation scores from BCI mortality survey subplots (Ramos et al., 2022; Schnitzer et al., 2021) |
| Fieldwork 2025 | Climbing mode and infestation pathway recorded in the field |

The code used to select the 33 TLS subplot locations is available at:  
https://github.com/PanamaForestGEO/TLSchoiceBCI

The following data are not included in this repository but are available upon reasonable request:

- 581 individually segmented tree point clouds, each labelled with ForestGEO tag number and species identity
- Individual subplot-level TLS point clouds (co-registered)
- Other raw data used in this analysis 
