context("DAPC")

test_that("test DAPC",{

  # inputs
  input <- dapc_items <- list()
  input$dapc_ploidy <- 2
  input$dapc_kmax <- 5
  input$dosage_file1$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")
  input$color_choice <- "YlOrRd"

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
  dapc_items$grp <- get_k$grp
  dapc_items$bestK <- get_k$bestK
  dapc_items$BIC <- get_k$BIC

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
  dapc_items$assignments <- clusters$Q
  dapc_items$dapc <- clusters$dapc


  ## Plots
  BIC <- dapc_items$BIC
  selected_K <- as.numeric(dapc_items$bestK)
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


  palette <- brewer.pal(as.numeric(input$dapc_k), input$color_choice)
  my_palette <- colorRampPalette(palette)(as.numeric(input$dapc_k))

  sc1 <- scatter.dapc(dapc_items$dapc,
                      bg = "white", solid = 1, cex = 1, # cex circle size
                      col = my_palette,
                      pch = 20, # shapes
                      cstar = 1, # 0 or 1, arrows from center of cluster
                      cell = 2, # size of elipse
                      scree.da = T, # plot da
                      scree.pca = T, # plot pca
                      posi.da = "topright",
                      posi.pca="bottomright",
                      mstree = F, # lines connecting clusters
                      lwd = 1, lty = 2,
                      leg = F, clab = 1) # legend and label of legend clusters. clab 0 or 1
})
