<!-- badges: start -->
[![R-CMD-check](https://github.com/Breeding-Insight/BIGapp/workflows/R-CMD-check/badge.svg)](https://github.com/Breeding-Insight/BIGapp/actions)
[![Development](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
  <!-- badges: end -->

# (B)reeding (I)nsight (G)enomics app <img src="https://github.com/user-attachments/assets/136c5ec4-5093-4129-a41b-233945e54198" align="right" width="250"/>

The BIGapp is a user-friendly tool for processing low to mid-density genotyping data for diploid and polyploid species. This R shiny app provides a web-based user friendly way for users to analyze genomic data without needing to use command-line tools. Additional analysis will be added, with the initial focus on a core set of features for supporting breeding decisions.

### Supported Analyses

Initial supported analyses will include the mature genomics/bioinformatics pipelines developed within Breeding Insight, with additional analyses continuing to be added.

Supported:
- Genotype processing
  - Dosage call from read counts
  - SNP filtering
  - Sample filtering
- Summary metrics
  - SNP Polymorphism Information Content
  - SNP Minor Allele Frequency
  - Sample Observed Heterozygosity
- Population Structure
  - PCA
  - DAPC
- GWAS
  - GWASpoly
- GS
  - Estimate Model Prediction Accuracy
  - Predict Phenotype Values and EBVs for Samples

### Running the BIG app

Tutorial available: 
https://scribehow.com/page/BIGapp_Tutorials__FdLsY9ZxQsi6kgT9p-U2Zg

Online preview:
https://big-demo.shinyapps.io/bigapp/

**Local computer**
1. Install R
2. Open Terminal (on mac)
3. To install and run development version of package:
(in terminal)
```
install.packages("devtools") #If not already installed
devtools::install_github("Breeding-Insight/BIGapp")
BIGapp::run_app()
```
4. View shiny app in browser

**Online (in progress)**

## Third-party software

The BIG app relies on both custom scripts and previously developed R packages cited below:

* [R](): version 4.2.2

#### R packages

* Shiny tools: [shiny](https://cran.r-project.org/web/packages/shiny/index.html), [shinyWidgets](https://cran.r-project.org/web/packages/shinyWidgets/index.html), [shinyalert](https://cran.r-project.org/web/packages/shinyalert/index.html), [shinyjs](https://cran.r-project.org/web/packages/shinyjs/index.html), [shinydisconnect](https://cran.r-project.org/web/packages/shinydisconnect/index.html), [shinycssloaders](https://cran.r-project.org/web/packages/shinycssloaders/index.html), [bs4Dash](https://cran.r-project.org/web/packages/bs4Dash/index.html),  [DT](https://cran.r-project.org/web/packages/DT/index.html), [config](https://cran.r-project.org/web/packages/config/index.html)

* Genetic analysis: [updog](https://cran.r-project.org/web/packages/updog/index.html), [GWASpoly](https://github.com/jendelman/GWASpoly), [AGHmatrix](https://cran.r-project.org/web/packages/AGHmatrix/index.html), [rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/index.html), [BIGr](https://github.com/Breeding-Insight/BIGr), [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html), [vcfR](https://cran.r-project.org/web/packages/vcfR/index.html)

* Data manipulation optimization: [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html), [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html), [purrr](https://cran.r-project.org/web/packages/purrr/index.html), [stringr](https://cran.r-project.org/web/packages/stringr/index.html), [future](https://cran.r-project.org/web/packages/future/index.html), [tibble](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html)

* Statistical analysis: [factoextra](https://cran.r-project.org/web/packages/factoextra/index.html), [MASS](https://cran.r-project.org/web/packages/MASS/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [matrixcalc](https://cran.r-project.org/web/packages/matrixcalc/index.html)

* Generate pretty graphics: [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html), [scales](https://cran.r-project.org/web/packages/scales/index.html), [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html), [plotly](https://cran.r-project.org/web/packages/plotly/index.html)
    

## Funding Sources
Breeding Insight is funded by USDA through Cornell University.
