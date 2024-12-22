context("diversity")

test_that("test diversity",{
  #Input variables (need to add support for VCF file)
  input <- diversity_items <- list()
  input$diversity_ploidy <- 2
  input$diversity_file$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")
  input$zero_value  <-  "Reference Allele Counts"
  input$hist_bins <- 20

  ploidy <- as.numeric(input$diversity_ploidy)
  geno <- input$diversity_file$datapath

  #Import genotype information if in VCF format
  vcf <- read.vcfR(geno)

  #Get items in FORMAT column
  info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT

  # Apply the function to the first INFO string
  info_ids <- extract_info_ids(info[1])

  #Get the genotype values if the updog dosage calls are present
  if ("UD" %in% info_ids) {
    geno_mat <- extract.gt(vcf, element = "UD")
    class(geno_mat) <- "numeric"
    rm(vcf) #Remove vcf
  }else{
    #Extract GT and convert to numeric calls
    geno_mat <- extract.gt(vcf, element = "GT")
    geno_mat <- apply(geno_mat, 2, convert_to_dosage)
    rm(vcf) #Remove VCF
  }

  #Convert genotypes to alternate counts if they are the reference allele counts
  #Importantly, the dosage plot is based on the input format NOT the converted genotypes
  is_reference <- (input$zero_value == "Reference Allele Counts")

  ######Get MAF plot (Need to remember that the VCF genotypes are likely set as 0 = homozygous reference, where the dosage report is 0 = homozygous alternate)

  #Status
  # Calculate percentages for both genotype matrices
  percentages1 <- calculate_percentages(geno_mat, ploidy)
  # Combine the data matrices into a single data frame
  percentages1_df <- as.data.frame(t(percentages1))
  percentages1_df$Data <- "Dosages"
  # Assuming my_data is your dataframe
  melted_data <- percentages1_df %>%
    pivot_longer(cols = -(Data),names_to = "Dosage", values_to = "Percentage")

  diversity_items$dosage_df <- melted_data

  #Convert the genotype calls prior to het,af, and maf calculation
  geno_mat <- data.frame(convert_genotype_counts(df = geno_mat, ploidy = ploidy, is_reference),
                         check.names = FALSE)

  # Calculating heterozygosity for a tetraploid organism
  diversity_items$het_df <- calculate_heterozygosity(geno_mat, ploidy = ploidy)

  diversity_items$maf_df <- calculateMAF(geno_mat, ploidy = ploidy)

  # Plots
  box <- ggplot(diversity_items$dosage_df, aes(x=Dosage, y=Percentage, fill=Data)) +
    #geom_point(aes(color = Data), position = position_dodge(width = 0.8), width = 0.2, alpha = 0.5) +  # Add jittered points
    geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.9) +
    labs(x = "\nDosage", y = "Percentage\n", title = "Genotype Distribution by Sample") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    )

  hist(diversity_items$het_df$Ho, breaks = as.numeric(input$hist_bins), col = "tan3", border = "black", xlim= c(0,1),
       xlab = "Observed Heterozygosity",
       ylab = "Number of Samples",
       main = "Sample Observed Heterozygosity")

  axis(1, at = seq(0, 1, by = 0.1), labels = TRUE)

  hist(diversity_items$maf_df$AF, breaks = as.numeric(input$hist_bins), col = "grey", border = "black", xlab = "Alternate Allele Frequency",
       ylab = "Frequency", main = "Alternate Allele Frequency Distribution")

  hist(diversity_items$maf_df$MAF, breaks = as.numeric(input$hist_bins), col = "grey", border = "black", xlab = "Minor Allele Frequency (MAF)",
       ylab = "Frequency", main = "Minor Allele Frequency Distribution")
})
