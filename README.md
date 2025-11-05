# Biomarkers and Pain Analysis


## Project Overview

This project investigates the relationship between inflammatory biomarkers and pain outcomes in patients recruited from two Scandinavian university hospitals.  
The analysis addresses two main research questions:

1. **Task 1 — Hypothesis Testing:** Do baseline inflammatory biomarker levels differ between male and female patients?  
2. **Task 2 — Regression Modelling:** Can baseline biomarkers and covariates predict 12-month pain outcomes?

The study uses **R** for all stages of analysis, including data preparation, hypothesis testing, regression modelling, and visualization.



## Project Structure
```
IPS-assessment/
│
├── README.md   # Project overview
│
├── main.r   # Main script for analysis
│
├── data/ 
│ ├── biomarkers.xlsx   # biomarker measurements
│ └── covariates.xlsx   # Demographics and clinical covariates
│
├── figs/   # Generated figures
  ├── dignostic_full.pdf
  ├── diagnostic_reduced.pdf
  ├── model_selection_detailed.pdf
  ├── model_selection.pdf
  ├── prediction_full.pdf
  └── prediction_reduced.pdf
```


## Requirements

The analysis was conducted in **R (version ≥ 4.3)**.  
The following packages are required:

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "janitor", "ggplot2",
  "broom", "car", "caret"
))
```


## Running the Analysis

1. **Open RStudio** (or R console) and set the working directory to this project folder.  

2. **Run the analysis scripts sequentially in the following order:**

```r
source("main.R")
```

3. **Generated figures will be saved in the `figs/` folder.**
- `dignostic_full.pdf`: Diagnostic plots for all biomarkers
- `diagnostic_reduced.pdf`: Diagnostic plots for selected biomarkers
- `model_selection.pdf`: Model selection plots for selected biomarkers
- `model_selection_detailed.pdf`: Model selection plots for all biomarkers (with detailed plots)
- `prediction_full.pdf`: Prediction plots for all biomarkers
- `prediction_reduced.pdf`: Prediction plots for selected biomarkers


## Key Results

- **Task 1:** No statistically significant sex differences in biomarker levels were found after Bonferroni correction, although **CSF-1** and **CXCL1** showed suggestive patterns that may warrant further investigation.

- **Task 2:** **IL-6**, **OPG**, and **baseline VAS** were statistically significant predictors of 12-month pain, but overall model performance remained modest (**Adjusted R² ≈ 0.27–0.33**), indicating limited predictive utility for clinical application.
