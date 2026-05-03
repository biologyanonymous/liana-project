# BCI Liana × Tree Architecture — Analysis Dataframes

**Project:** Tree Architectural Feedbacks with Liana Canopy Infestation  

These six CSVs are the final, analysis-ready datasets for the BCI thesis project. They are produced by `save_final_dfs.R` (run once after `final_figures_v13.R`) and are read directly by `bci_liana_architecture.qmd`. No raw-data joining or cleaning is needed in the analysis script.

---

## Files

### `dat_h1.csv` — H1 structural legacies dataset
**n = 251 trees | 67 species | 24 subplots**

Trees from the BCI 50-ha plot with TLS crown metrics (2019) and infestation histories reconstructed from annual dendrometer records (2011–2019). Filtered to three mutually exclusive infestation-history categories; recently infested and uncertain trees excluded.

| Column | Description |
|--------|-------------|
| `tag6` | 6-digit ForestGEO stem tag |
| `species` | Species name (factor) |
| `sbpltnm` | 40×40 m subplot ID (factor) |
| `H1_category` | `Never infested` / `Lost lianas` / `Persistently infested` |
| `dbh_2019` | DBH in mm (raw dendrometer measurement) |
| `dbh_cm` | DBH in cm (`dbh_2019 / 10`) |
| `height` | Total tree height in m (TLS 2019) |
| `crown_depth` | Crown depth in m (TLS) |
| `pa` | Crown projected area in m² (TLS) |
| `crown_vol` | Crown volume in m³ (TLS) |
| `crown_base_h` | Height of crown base above ground (`height − crown_depth`) |
| `slenderness` | Slenderness index (`height / dbh_cm`) |
| `cbh_ratio` | Crown base height ratio (`crown_base_h / height`) |
| `log_dbh` | `log(dbh_cm)` |
| `log_height` | `log(height)` |
| `log_depth` | `log(crown_depth)` |
| `log_area` | `log(pa)` |
| `log_vol` | `log(crown_vol)` |
| `log_slender` | `log(slenderness)` |
| `log_cbh` | `log(max(cbh_ratio, 0.001))` — floor prevents log(0) |

**Used by:** H1 Bayesian LMMs (brms), slenderness/CBH violin plots, 4-panel allometric figure, species RE scatter plot.

---

### `clean_pai.csv` — H1 plant area index dataset
**n = 575 trees | 91 species | 24 subplots**

Voxel-based plant area index (PAI) estimated from TLS point clouds. Liana material extracted using alpha-shape envelopes; PAI calculated separately for host-tree points and combined host+liana points. Joined with subplot IDs from the dendrometer records.

| Column | Description |
|--------|-------------|
| `tag` | Tree tag (numeric) |
| `Species` | Species name (factor) |
| `sbpltnm` | Subplot ID (factor) |
| `DBH` | DBH in mm |
| `dbh_cm` | DBH in cm |
| `log_dbh` | `log(dbh_cm)` |
| `Lianas` | Liana cover score (0–100%) from 2019 field survey |
| `infest_cat` | Liana load category: `None` / `Light` (1–50%) / `Heavy` (51–100%) |
| `infest_bin` | Binary: `None` / `Infested` |
| `PAI_tree_only` | Host PAI (TLS, tree points only) |
| `PAI_combined` | Combined PAI (tree + liana points) |

**Used by:** H1 PAI Bayesian LMM (brms), all PAI bar charts including top-12 species appendix figure.

---

### `dat_h2.csv` — H2 pre-infestation predictor dataset
**n = 181 trees | 61 species | 24 subplots**

Subset of trees uninfested in 2019, classified by whether they subsequently gained liana infestation by 2025 (`gained infestation`) or remained uninfested (`remained uninfested`). Crown metrics are TLS 2019 measurements — i.e. the baseline architecture *before* the outcome is determined. Column structure is identical to `dat_h1.csv` plus `H2_category`.

| Column | Description |
|--------|-------------|
| `H2_category` | `remained uninfested` / `gained infestation` |
| *(all others)* | Same as `dat_h1.csv` |

**Used by:** H2 Bayesian LMMs (brms), H2 4-panel figure, H2 height-by-species violin plot.

---

### `dat_mode.csv` — H2b climbing mode × crown architecture dataset

Fieldwork data collected in 2025 for the 6 most abundant focal species (≥10 individuals per species), joined with TLS crown metrics. One row per tree. Dominant climbing mode and class already summarised from potentially multiple observed modes per tree.

| Column | Description |
|--------|-------------|
| `tag6` | 6-digit ForestGEO stem tag |
| `species` | Species name (factor) |
| `infestation_type` | `Direct` / `Lateral` / `Both` (how lianas entered the crown) |
| `dominant_mode` | Most common climbing mode observed on that tree |
| `dominant_class` | `Active` (twining/tendrils/prehensile) or `Passive` (scrambling/hooks/adhesive) |
| `crown_depth`, `pa`, `crown_vol`, `height` | TLS crown metrics (m, m², m³, m) |
| `dbh_cm`, `log_dbh` | Size covariate |
| `log_depth`, `log_area`, `log_vol`, `log_height` | Log-transformed metrics for OLS models |

`dat_mode_ap` is derived from this file in the analysis script by filtering to rows where `dominant_class` is not NA.

**Used by:** H2b OLS models (lm), Active vs Passive crown architecture boxplots.

---

### `pathway_dat.csv` — infestation pathway × H2 outcome dataset

Fieldwork records with infestation pathway type joined with H2 infestation outcome. Used for the chi-square test of whether pathway type (Direct/Lateral/Both) is associated with whether a tree subsequently gains or retains infestation.

| Column | Description |
|--------|-------------|
| `tag6` | 6-digit ForestGEO stem tag |
| `H2_category` | Infestation outcome (factor: `gained infestation` / `remained uninfested`) |
| `infestation_type` | `Direct` / `Lateral` / `Both` |
| `inf_type_bin` | Binary: `Lateral` vs `Non-lateral` |
| `lateral` | Integer (0/1): whether liana entry was lateral |

**Used by:** §7 chi-square test, infestation-type-by-species proportional bar chart.

---

### `modes_clean.csv` — climbing mode records for species bar charts

Long-format fieldwork records: one row per tree × observed liana climbing mode (a tree with multiple liana species can appear multiple times). Used for the species-level proportional bar charts of climbing mode and class.

| Column | Description |
|--------|-------------|
| `tag6` | 6-digit ForestGEO stem tag |
| `species` | Species name |
| `mode` | Climbing mode (e.g. `twining`, `scrambling`, `hooks`, `tendrils`) |
| `infestation_type` | `Direct` / `Lateral` / `Both` |
| `climb_class` | `Active` or `Passive` (derived from `mode`) |

**Used by:** §8.3 climbing mode and class proportional bar charts (top 6 focal species only).

---

## Reproducing the analysis

```
~/mbiol-project/
├── data/
│   ├── dat_h1.csv
│   ├── clean_pai.csv
│   ├── dat_h2.csv
│   ├── dat_h2_2011.csv
│   ├── dat_mode.csv
│   ├── pathway_dat.csv
│   └── modes_clean.csv
├── models/               ← pre-fitted brms .rds files (load instantly; delete to refit)
│   ├── m_h1.rds
│   ├── m_h2_lmm.rds
│   └── m_pai_main.rds
├── figures/
├── tables/
├── bci_liana_architecture.qmd   ← main analysis document
├── save_final_dfs.R             ← script that produced these files
└── README.md                    ← this file
```

Open `bci_liana_architecture.qmd` in RStudio and click **Render**. All figures are written to `figures/` and tables to `tables/` automatically.

To refit a Bayesian model from scratch, delete the corresponding `.rds` file from `models/` and re-render.

---

## Data provenance

| Source | Description |
|--------|-------------|
| TLS 2019 | Point clouds acquired across 24 BCI 50-ha plot subplots; processed with ITSMe (Terryn et al., 2023) and RayCloudTools (Lowe & Stepanas, 2021) |
| Dendrometer records 2011–2025 | Annual liana infestation scores (0–4 scale) from BCI mortality survey subplots (Ramos et al., 2022; Schnitzer et al., 2021) |
| Fieldwork 2025 | Climbing mode and infestation pathway recorded for 6 focal species (≥10 individuals each) |

Raw data and processing scripts are held separately and are available on request.
