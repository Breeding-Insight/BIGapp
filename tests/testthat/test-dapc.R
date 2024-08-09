context("DAPC")

test_that("test DAPC",{

  # inputs
  input <- list()
  input$dapc_ploidy <- 2
  input$dapc_kmax <- 5
  input$dosage_file1$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")

  input$plot_BICX <- TRUE

  # Step 1
  ploidy <- as.numeric(input$dapc_ploidy)
  maxK <- as.numeric(input$dapc_kmax)
  geno <- input$dosage_file1$datapath

  ##Add in VCF with the vcfR package (input VCF, then convert to genlight using vcf2genlight function)

  #Import genotype information if in VCF format
  vcf <- read.vcfR(geno)

  #Get items in FORMAT column
  info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT

  # Apply the function to the first INFO string
  info_ids <- extract_info_ids(info[1])

  #Get the genotype values if the updog dosage calls are present
  if ("UD" %in% info_ids) {
    genotypeMatrix <- extract.gt(vcf, element = "UD")
    class(genotypeMatrix) <- "numeric"
    rm(vcf) #Remove vcf
  }else{
    #Extract GT and convert to numeric calls
    genotypeMatrix <- extract.gt(vcf, element = "GT")
    genotypeMatrix <- apply(genotypeMatrix, 2, convert_to_dosage)
    rm(vcf) #Remove VCF
  }

  #Perform analysis
  get_k <- findK(genotypeMatrix, maxK, ploidy)

  #Assign results to reactive values
  get_k # Containg bestK, grp, BIC
  dapc_items1 <- get_k

  # Step 2
  # Inputs
  input <- list()
  input$dapc_ploidy <- 2
  input$dapc_k <- 5
  input$dosage_file2$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")

  geno <- input$dosage_file2$datapath
  ploidy <- as.numeric(input$dapc_ploidy)
  selected_K <- as.numeric(input$dapc_k)

  #Import genotype information if in VCF format
  vcf <- read.vcfR(geno)

  #Get items in FORMAT column
  info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT

  # Apply the function to the first INFO string
  info_ids <- extract_info_ids(info[1])

  #Get the genotype values if the updog dosage calls are present
  if ("UD" %in% info_ids) {
    genotypeMatrix <- extract.gt(vcf, element = "UD")
    class(genotypeMatrix) <- "numeric"
    rm(vcf) #Remove vcf
  }else{
    #Extract GT and convert to numeric calls
    genotypeMatrix <- extract.gt(vcf, element = "GT")
    genotypeMatrix <- apply(genotypeMatrix, 2, convert_to_dosage)
    rm(vcf) #Remove VCF
  }

  #Perform analysis
  clusters <- performDAPC(genotypeMatrix, selected_K, ploidy)

  #Assign results to reactive value
  clusters # contains Q and dapc
  dapc_items2 <- clusters


  ## Plots
  BIC <- dapc_items1$BIC
  selected_K <- as.numeric(dapc_items1$bestK)
  plot(BIC, type = "o", xaxt = 'n')
  axis(1, at = seq(1, nrow(BIC), 1), labels = TRUE)

  if (input$plot_BICX) {
    p <- plot(BIC, type = "o", xaxt = 'n')
    axis(1, at = seq(1, nrow(BIC), 1), labels = TRUE)
    points(selected_K, BIC[selected_K,2], pch = "x", col = "red", cex = 2)
  } else {
    plot(BIC, type = "o", xaxt = 'n')
    axis(1, at = seq(1, nrow(BIC), 1), labels = TRUE)
  }

})
