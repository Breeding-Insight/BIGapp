# (B)reeding (I)nsight (G)enomics app

Currently, Breeding Insight provides bioinformatic processing support for our external collaborators. This R shiny app will provide a web-based user friendly way for our internal and external collaborators to analyze genomic data without needing to use command-line tools.

### Supported Analyses
Initial supported analyses will include the mature genomics/bioinformatics pipelines developed within Breeding Insight, with additional analyses continuing to be added.

Supported:
- Genotype processing
  - Dosage call from read counts
  - SNP filtering
  - Sample filtering
- Summary metrics
  - SNP Allele Frequency
  - SNP Minor Allele Frequency
  - Sample Observed Heterozygosity
- Population Structure
  - PCA
  - DAPC
- GWAS
- GS

### Running the BIG app

**Local computer**

1. Install R
2. Download this folder from GitHub
3. Open Terminal (mac)
4. Initialize R in Terminal and enter the below commands
```
install.packages("shiny") #This is only needed for the first time running shiny
library("shiny") #Starting the shiny package
setwd("path_to_shiny_app/Genomics_Shiny_App") #Direct R to the downloaded app folder
runApp("BIG_app") #Start the app
```
8. View shiny app in browser

**Online (in progress)**

## References
The BIG app relies on both custom scripts and previously developed R packages cited below:

R: version 4.2.2

required_cran_packages <- c("updog", "ggplot2","devtools","GWASpoly","SNPRelate",
                       "adegenet", "future", "scales", "AGHmatrix", "stats", 
                       "factoextra", "readxl", "ggrepel", "dplyr", "shiny",
                       "shinydashboard","randomcoloR","plotly", "DT","RColorBrewer",
                       "dichromat", "bs4Dash", "shinyWidgets","data.table",
                       "matrixcalc","Matrix", "shinyalert","rrBLUP", "tidyverse",
                       "foreach", "doParallel","VariantAnnotation", "vcfR")

required_bio_packages <- c("SNPRelate","VariantAnnotation")

Dev_tools_packages <- c("GWASpoly")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    

## Funding Sources
Breeding Insight is funded by USDA through Cornell University.
