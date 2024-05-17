# (B)reeding (I)nsight (G)enomics app

Currently, Breeding Insight provides bioinformatic processing support for our external collaborators. This app will a provide web-based user friendly way for our internal and external collaborators to analyze genomic data without needing to use command-line tools.

### Supported Analyses
Initial supported analyses will include the mature genomics/bioinformatics pipelines developed within Breeding Insight, with additional analyses continuing to be added.

Supported:
1. Dosage calling from MADC file
2. SNP filtering
3. Population Structure
4. Genomic Diversity
5. GWAS

In-progress:
1. Genomic Prediction/Selection

### Running the BIG app online (in progress)



### Running the BIG app on a local computer
To use:
1. Install.packages("shiny")
2. library("shiny")
3. setwd("/Genomics_Shiny_App")
4. runApp("BIG_app")
5. View shiny app in browser

## References
The BIG app relies on both custom scripts and previously developed R packages cited below:

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
    
**jhg**

## Funding Sources
