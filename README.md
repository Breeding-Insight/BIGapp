<!-- badges: start -->
[![R-CMD-check](https://github.com/Breeding-Insight/BIGapp/workflows/R-CMD-check/badge.svg)](https://github.com/Breeding-Insight/BIGapp/actions)
![GitHub Release](https://img.shields.io/github/v/release/Breeding-Insight/BIGapp)
[![Development Status](https://img.shields.io/badge/development-active-blue.svg)](https://img.shields.io/badge/development-active-blue.svg)
![GitHub License](https://img.shields.io/github/license/Breeding-Insight/BIGapp)



<!-- badges: end -->

<div align="center">

# (B)reeding (I)nsight (G)enomics app (BIGapp)

<img src="https://github.com/user-attachments/assets/136c5ec4-5093-4129-a41b-233945e54198"  width="250"/>

</div>

BIGapp is a user-friendly web application built with R and Shiny, designed to simplify the processing of low to mid-density genotyping data for both diploid and polyploid species. It provides a powerful and intuitive interface for researchers and breeders to analyze genomic data without requiring command-line expertise.

## Key Features

- **Web-Based Interface:** Access BIGapp through your web browser, eliminating the need for using command-line inputs to perform genomic analysis.
- **Genotype Processing:**
    -  Call genotypes from read counts.
    -  Filter SNPs based on various criteria.
    -  Filter samples to ensure data quality.
- **Summary Statistics:**
    -  Calculate SNP Polymorphism Information Content (PIC).
    -  Determine SNP Minor Allele Frequency (MAF).
    -  Compute Sample Observed Heterozygosity.
- **Population Structure Analysis:**
    -  Perform Principal Component Analysis (PCA).
    -  Conduct Discriminant Analysis of Principal Components (DAPC).
- **Genome-Wide Association Studies (GWAS):**
    -  Utilize GWASpoly for robust association mapping.
- **Genomic Selection (GS):**
    -  Estimate model prediction accuracy.
    -  Predict phenotypic values and Estimated Breeding Values (EBVs) for your samples.
- **Expanding Functionality:** BIGapp is actively developed, with new analyses and features continuously being added.

## User Interface

<p align="center">
  <img src="https://github.com/user-attachments/assets/9a6984df-8116-403c-85c1-ba9600623940" alt="BIGapp Screenshot" width="800"/>
  <br>
  <em>BIGapp's intuitive interface makes genomic data analysis accessible to everyone.</em>
</p>

## Getting Started

### Tutorials
New to BIGapp? Check out our comprehensive tutorial to guide you through the process: [BIGapp Tutorials](https://scribehow.com/page/BIGapp_Tutorials__FdLsY9ZxQsi6kgT9p-U2Zg)

### Online Preview

Try out a live demo of BIGapp here: [BIGapp Demo](https://big-demo.shinyapps.io/bigapp-main/)

### Local Installation

1. **Install R:** Download and install the latest version of R from [CRAN](https://cran.r-project.org/).
2. **Open Terminal (macOS/Linux) or R Console (Windows).**
3. **Installation:**
    ```R
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        install.packages("remotes")
    
    BiocManager::install("Breeding-Insight/BIGapp", dependencies = TRUE)
    ```
4. **Starting BIGapp:**
    ```R
    BIGapp::run_app()
    ```
5. **Access in Browser:** The BIGapp interface will open in your default web browser.

### Online Deployment (Coming Soon)

BIGapp will be deployed on USDA SciNet for convenient online access. Stay tuned for updates!

## Dependencies

BIGapp leverages a powerful suite of R packages:

### Core R Packages

-   **R (>= 4.4.0)**

### Shiny Framework

-   [shiny](https://cran.r-project.org/web/packages/shiny/index.html): Web application framework.
-   [shinyWidgets](https://cran.r-project.org/web/packages/shinyWidgets/index.html): Custom input widgets.
-   [shinyalert](https://cran.r-project.org/web/packages/shinyalert/index.html): Create elegant pop-up messages.
-   [shinyjs](https://cran.r-project.org/web/packages/shinyjs/index.html): Enhance Shiny apps with JavaScript actions.
-   [shinydisconnect](https://cran.r-project.org/web/packages/shinydisconnect/index.html): Handle user disconnections gracefully.
-   [shinycssloaders](https://cran.r-project.org/web/packages/shinycssloaders/index.html): Add CSS loaders for visual feedback.
-   [bs4Dash](https://cran.r-project.org/web/packages/bs4Dash/index.html): Bootstrap 4 dashboard components.
-   [DT](https://cran.r-project.org/web/packages/DT/index.html): Display data tables with interactive features.
-   [config](https://cran.r-project.org/web/packages/config/index.html): Manage environment-specific configurations.

### Genetic Analysis

-   [updog](https://cran.r-project.org/web/packages/updog/index.html): Genotype polyploid individuals.
-   [GWASpoly](https://github.com/jendelman/GWASpoly): Conduct GWAS in polyploids.
-   [AGHmatrix](https://cran.r-project.org/web/packages/AGHmatrix/index.html): Compute genomic relationship matrices.
-   [rrBLUP](https://cran.r-project.org/web/packages/rrBLUP/index.html): Perform genomic prediction.
-   [BIGr](https://github.com/Breeding-Insight/BIGr): Breeding Insight's core genomic analysis functions.
-   [adegenet](https://cran.r-project.org/web/packages/adegenet/index.html): Explore and analyze genetic data.
-   [vcfR](https://cran.r-project.org/web/packages/vcfR/index.html): Manipulate and analyze VCF files.

### Data Manipulation

-   [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html): Data manipulation tools.
-   [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html): Tidy your data.
-   [purrr](https://cran.r-project.org/web/packages/purrr/index.html): Functional programming toolkit.
-   [stringr](https://cran.r-project.org/web/packages/stringr/index.html): String manipulation.
-   [future](https://cran.r-project.org/web/packages/future/index.html): Unified parallel and distributed processing.
-   [tibble](https://cran.r-project.org/web/packages/tibble/vignettes/tibble.html): Modern data frame alternative.

### Statistical Analysis

-   [factoextra](https://cran.r-project.org/web/packages/factoextra/index.html): Extract and visualize multivariate analyses.
-   [MASS](https://cran.r-project.org/web/packages/MASS/index.html): Statistical functions and datasets.
-   [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html): Sparse and dense matrix operations.
-   [matrixcalc](https://cran.r-project.org/web/packages/matrixcalc/index.html): Matrix calculus functions.

### Visualization

-   [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html): Create elegant data visualizations.
-   [scales](https://cran.r-project.org/web/packages/scales/index.html): Graphical scaling methods.
-   [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html): Color palettes for thematic maps.
-   [plotly](https://cran.r-project.org/web/packages/plotly/index.html): Create interactive web graphics.

## Funding

BIGapp development is supported by [Breeding Insight](https://www.breedinginsight.org/), a USDA-funded initiative based at Cornell University.

## Citation

If you use BIGapp in your research, please cite:

Sandercock A.M., Peel M.D., Tanigut C.H., Chinchilla-Vargas J., Chen S., Sapkota M., Lin M., Zhao D., Ackerman A.J., Basnet B.R., Beil C.T., Sheehan M.J. (2025). BIGapp: A User-Friendly Genomic Tool Kit Identified Quantitative Trait Loci for Creeping Rootedness in Alfalfa (_Medicago sativa_ L.)., The Plant Genome. doi:https://doi.org/10.1002/tpg2.70067
