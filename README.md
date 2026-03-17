# NHIS 2024 — GLP-1 Receptor Agonist Use Among Cancer Survivors With Diabetes

**Study:** Disparities in GLP-1 Receptor Agonist Use Among Adults With Diabetes and Cancer History: Findings From the 2024 National Health Interview Survey

**Authors:** Olinto Linares-Perdomo, PhD, MS, Damon Klebe, PhD, Bismarck C Odei, MD  
**Affiliation:** Psycho-Oncology, Department of Radiation Oncology, Huntsman Cancer Institute, University of Utah School of Medicine, Salt Lake City, UT

---

## Overview

This repository contains the complete R analysis pipeline for a cross-sectional study using the 2024 National Health Interview Survey (NHIS-2024). The study evaluates GLP-1 receptor agonist (GLP-1 RA) use among U.S. adults with diagnosed diabetes, comparing cancer survivors to those without cancer history, and assessing effect modification by insurance type and race/ethnicity.

**Key findings:**
- Cancer survivors had approximately 30% lower adjusted odds of GLP-1 RA use (aOR 0.70; 95% CI 0.52–0.96)
- Adjusted gap: −6.0 percentage points (95% CI −6.6 to −5.4)
- Insurance type significantly modified the association (P for interaction = 0.045)

---

## Repository Structure

```
nhis2024-glp1-cancer-survivorship/
├── README.md
├── LICENSE
├── sessionInfo.txt          # R and package versions used
├── .gitignore
├── data/
│   └── README.md            # Data download instructions (data not included)
├── R/
│   ├── Block_00_data_loading.R
│   ├── Block_01_variable_construction.R
│   ├── Block_02_survey_design.R
│   ├── Block_03_table1_descriptives.R
│   ├── Block_04_table2_primary_analysis.R
│   ├── Block_05_table3_effect_modification.R
│   ├── Block_06_etable1_regression.R
│   ├── Block_07_etable2_sensitivity.R
│   ├── Block_08_etable3_model_c.R
│   ├── Block_09_figure1.R
│   └── Block_10_supplemental_race_table.R
|   └── run_all.R
└── output/
    └── README.md            # Description of pipeline outputs
```

---

## Data Provenance

This pipeline uses the **NHIS 2024 Sample Adult file**, a publicly available dataset produced by the National Center for Health Statistics (NCHS).

**Download instructions:**

1. Go to: https://www.cdc.gov/nchs/nhis/2024nhis.htm
2. Download the **Sample Adult** CSV file (`adult24.csv`)
3. Place `adult24.csv` in the `data/` folder at the project root

The data are free to download and do not require registration. Full documentation and variable codebooks are available at the same URL.

**Note:** The data file is not included in this repository in accordance with NCHS terms of use. The pipeline will check for the file at startup and provide a download reminder if it is not found.

---

## Requirements

### R version
R 4.4.0 or later (see `sessionInfo.txt` for exact version used)

### Required packages

```r
install.packages(c(
  "tidyverse",   # data manipulation
  "survey",      # complex survey analysis (version 4.2-1 or later)
  "flextable",   # table formatting
  "officer",     # Word document export
  "ggplot2",     # figures
  "here"         # project-relative file paths
))
```

See `sessionInfo.txt` for exact package versions used in the published analysis.

---

## How to Run

Blocks must be run **in order** (00 → 10). Each block depends on objects created by the preceding block.

```r
# Option 1: Run all blocks sequentially
source("R/Block_00_data_loading.R")
source("R/Block_01_variable_construction.R")
source("R/Block_02_survey_design.R")
source("R/Block_03_table1_descriptives.R")
source("R/Block_04_table2_primary_analysis.R")
source("R/Block_05_table3_effect_modification.R")
source("R/Block_06_etable1_regression.R")
source("R/Block_07_etable2_sensitivity.R")
source("R/Block_08_etable3_model_c.R")
source("R/Block_09_figure1.R")
source("R/Block_10_supplemental_race_table.R")

For convenience, the full pipeline can be executed in a single step:
# Option 2: Run full pipeline (recommended)
source("R/run_all.R")
```


Alternatively, open the project in RStudio and run blocks interactively.

**Recommended:** Use the `here` package for project-relative paths (already implemented). Open the `.Rproj` file or set your working directory to the project root before running.

---

## Block Descriptions

| Block | File | Output |
|-------|------|--------|
| 00 | `Block_00_data_loading.R` | `final_df` — analytic sample (n = 3,615) |
| 01 | `Block_01_variable_construction.R` | `analytic_df` — recoded analysis variables |
| 02 | `Block_02_survey_design.R` | `nhis_design` — svydesign object |
| 03 | `Block_03_table1_descriptives.R` | `table1_baseline_characteristics.docx` |
| 04 | `Block_04_table2_primary_analysis.R` | `table2_glp1_models.docx`; `model_a`, `model_b` |
| 05 | `Block_05_table3_effect_modification.R` | `table3_effect_modification.docx` |
| 06 | `Block_06_etable1_regression.R` | `eTable1_regression_coefficients.docx` |
| 07 | `Block_07_etable2_sensitivity.R` | `eTable2_sensitivity_analyses.docx` |
| 08 | `Block_08_etable3_model_c.R` | `eTable3_model_c_attenuation.docx`; `model_c` |
| 09 | `Block_09_figure1.R` | `Figure1_adjusted_probabilities.png/.pdf` |
| 10 | `Block_10_supplemental_race_table.R` | `GLP1_supplemental_race_table.docx` |

---

## Analytic Methods

### Study design
Cross-sectional analysis of the NHIS 2024 Sample Adult file, a nationally representative survey of the U.S. civilian noninstitutionalized population.

### Primary exposure
Cancer history: any self-reported clinician-diagnosed cancer, excluding nonmelanoma skin cancer only (consistent with NCI SEER conventions).

### Primary outcome
Current GLP-1 RA use: NHIS variable DIBGLP_A, new Emerging Content in the 2024 cycle.

### Statistical approach
- Survey-weighted logistic regression with Taylor series linearization
- Marginal standardization (G-computation) for adjusted predicted probabilities
- Interaction tests using Wald tests (regTermTest) at α = 0.10
- Cell suppression threshold: n ≥ 30 per cancer history group per stratum (NCHS guidelines)

### Pre-specified analytic adaptations
| Adaptation | Rationale |
|-----------|-----------|
| Insurance: 4 → 3 categories (uninsured excluded) | n = 6 uninsured cancer survivors; sparse for regression |
| BMI: 4 → 3 categories (<25 + 25–29.9 → <30) | Sparse cell at BMI <25 in cancer stratum |
| Race/ethnicity: 5 → 4 categories (NH Asian → "Other") | Sparse cells in cancer stratum |
| Outpatient visit excluded from regression | 98.5% prevalence; no discriminatory variation |


---

## Reproducibility

All results reported in the manuscript (tables, figures, and supplementary analyses) can be reproduced by running the R scripts in the `R/` directory sequentially.

The pipeline generates all analytic outputs in the `output/` directory.
---

## Reporting Guidelines

This study was reported following the STROBE guidelines for observational studies.

---

## Citation

If you use this pipeline, please cite the associated manuscript:

> Linares-Perdomo O, Klebe D, Odei BC. Disparities in GLP-1 Receptor Agonist Use Among Adults With Diabetes and Cancer History: Findings From the 2024 National Health Interview Survey. *[Journal]*. [Year]. doi:[DOI]

---

## License

This code is released under the MIT License. See `LICENSE` for details.

The NHIS 2024 data are in the public domain; see NCHS terms of use at:  
https://www.cdc.gov/nchs/nhis/about_nhis.htm

---

## Version

This repository corresponds to the analytic version used for manuscript submission.

Future updates may reflect revisions following peer review.

---

## Contact

Olinto Linares-Perdomo, PhD, MS  
Email: 
olinto_Linares@abl.med.utah.edu  
olinto.linares.perdomo@gmail.com

Psycho-Oncology, Department of Radiation Oncology  
Huntsman Cancer Institute, University of Utah School of Medicine  
1950 Circle of Hope Dr., Rm 1570, Salt Lake City, UT 84112-USA

