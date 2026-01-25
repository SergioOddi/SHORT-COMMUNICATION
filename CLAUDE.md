# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a scientific short communication project that combines data analysis in R with manuscript writing. The workflow involves:
1. Data analysis and visualization using R scripts
2. Scientific manuscript preparation with proper reference management
3. Integration of figures and results into the final document

## Project Structure

```
SHORT-COMMUNICATION/
├── DATASET/          # Raw and processed data files
├── SCRIPTS/          # R scripts for analysis and figure generation
├── FIGURES/          # High-resolution figures for publication
├── MANUSCRIPT/       # Scientific manuscript in R Markdown format
├── REFERENCES/       # Bibliography files (BibTeX format)
└── OUTPUT/           # Intermediate results, tables, statistics
```

## Workflow

### 1. Data Analysis
- R scripts should be numbered sequentially (e.g., `01_`, `02_`, `03_`)
- Scripts should be self-contained and well-commented
- Save intermediate results to OUTPUT/
- Generate figures and save to FIGURES/ in publication-ready format (PDF or high-res PNG/TIFF at 300 dpi)

### 2. Manuscript Preparation
- Use R Markdown (.Rmd) format for the manuscript
- Store in MANUSCRIPT/ directory
- Include YAML header with output format specifications
- Reference bibliography file: `../REFERENCES/bibliography.bib`

### 3. Reference Management
- Use BibTeX format (.bib file) in REFERENCES/
- Citation syntax in R Markdown: `[@citation_key]`
- Bibliography automatically generated at end of document

## Key Technologies

- **R**: Primary language for statistical analysis
- **R Markdown**: Manuscript format that integrates code, results, and narrative
- **knitr/rmarkdown**: For rendering Rmd to PDF/Word/HTML
- **BibTeX**: Reference management

## Common Commands

### Running R Scripts
```r
# In R console or RStudio
source("SCRIPTS/01_script_name.R")
```

### Rendering Manuscript
```r
# In R console
rmarkdown::render("MANUSCRIPT/manuscript.Rmd", output_format = "pdf_document")

# For Word output
rmarkdown::render("MANUSCRIPT/manuscript.Rmd", output_format = "word_document")
```

### Installing Required Packages
```r
# Check and install common packages
required_packages <- c("tidyverse", "ggplot2", "knitr", "rmarkdown")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
```

## R Markdown Manuscript Template Structure

```yaml
---
title: "Your Title"
author: "Author Name"
date: "`r Sys.Date()`"
output:
  pdf_document:
    fig_caption: yes
    keep_tex: yes
bibliography: ../REFERENCES/bibliography.bib
csl: journal-style.csl  # Optional: specific journal citation style
---
```

Main sections:
- Abstract
- Introduction
- Methods
- Results (with embedded figures)
- Discussion
- References (auto-generated)

## Figure Integration in R Markdown

```r
# Code chunk with figure generation
```{r figure-label, fig.cap="Caption text", fig.width=7, fig.height=5}
# R code to generate plot
```

## Citation Management

Add references to `REFERENCES/bibliography.bib`:
```bibtex
@article{author2023,
  title={Article Title},
  author={Author, A. and Author, B.},
  journal={Journal Name},
  year={2023},
  volume={10},
  pages={1-10}
}
```

Cite in text: `[@author2023]` or `@author2023` for inline citation.

## Language Note

Rispondere sempre in italiano quando si interagisce con l'utente.
