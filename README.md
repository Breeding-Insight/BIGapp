# (B)reeding (I)nsight (G)enomics app

Currently, Breeding Insight provides bioinformatic processing support for our external collaborators. This app will a provide web-based user friendly way for our internal and external collaborators to analyze genomic data without needing to use command-line tools.

### Supported Analyses
Initial supported analyses will include the mature genomics/bioinformatics pipelines developed within Breeding Insight, with additional analyses continuing to be added.

Supported:
- Dosage calling from MADC file
- SNP filtering
- Population Structure
- Genomic Diversity
- GWAS

In-progress:
- Genomic Prediction/Selection

### Running the BIG app

**Local computer**
1. Install R
2. Install.packages("shiny")
3. library("shiny")
4. setwd("/Genomics_Shiny_App")
5. runApp("BIG_app")
6. View shiny app in browser

**Online (in progress)**

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
