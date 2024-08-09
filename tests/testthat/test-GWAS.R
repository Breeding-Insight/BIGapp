context("GWAS")

test_that("test GWAS",{
  input <- list()
  input$cores <- 1
  input$phenotype_file$datapath <- system.file("iris_passport_file.csv", package = "BIGapp")
  input$trait_info <- "Sepal.Length"
  input$trait_info <- c("Sepal.Length", "Sepal.Width")
  input$fixed_info <- "Species"
  input$gwas_ploidy <- 2
  input$gwas_file$datapath <- system.file("iris_DArT_VCF.vcf.gz", package = "BIGapp")
  input$gwas_threshold <- "M.eff"
  input$model_select <- "additive"

  cores <- input$cores

  #Make subset phenotype file (need to develop alternative method that does not save a new phenotype file each time.)
  #I think I can subset the read.GWAS file pheno and fixed categories (data@pheno[,c("trait")]) and data@fixed = phenotype_file[,c("List of fixed traits")]
  phenotype_file <- read.csv(input$phenotype_file$datapath, header = TRUE, check.names = FALSE)

  ids <- colnames(phenotype_file)[1]
  traits <- input$trait_info
  fixed <- input$fixed_info
  included_var <- c(ids, traits, fixed)
  ploidy <- as.numeric(input$gwas_ploidy)

  # Check if traits are numerical
  n_traits <- as.matrix(phenotype_file[,traits])
  n_traits <- apply(n_traits, 2, function(x) all(is.na(as.numeric(x))))
  if(any(n_traits)) stop("The selected traits must be numerical.")

  phenotype_file <- phenotype_file[,included_var]

  # Create a temporary file for the selected phenotype data
  temp_pheno_file <- tempfile(fileext = ".csv")

  #Save new phenotype file with selected traits and fixed effects
  write.csv(phenotype_file, file = temp_pheno_file, row.names = FALSE)

  #Remove the phenotype_file from memory
  rm(phenotype_file)

  #Geno file path
  file_path <- input$gwas_file$datapath

  #Geno.file conversion if needed
  if (grepl("\\.csv$", file_path)) {
    data <- read.GWASpoly(ploidy= ploidy, pheno.file= temp_pheno_file, geno.file=input$gwas_file$datapath,
                          format="numeric", n.traits=length(traits), delim=",") #only need to change files here

  } else if (grepl("\\.vcf$", file_path) || grepl("\\.gz$", file_path)) {
    # Create a temporary file for the selected phenotype data
    temp_geno_file <- tempfile(fileext = ".csv")

    #Convert VCF file if submitted
    vcf <- read.vcfR(input$gwas_file$datapath)

    #Extract GT
    geno_mat <- extract.gt(vcf, element = "GT")
    geno_mat <- apply(geno_mat, 2, convert_to_dosage)
    class(geno_mat) <- "numeric"
    info <- data.frame(vcf@fix)
    gpoly_df <- cbind(info[,c("ID","CHROM","POS")], geno_mat)
    write.csv(gpoly_df, file = temp_geno_file, row.names = FALSE)

    data <- read.GWASpoly(ploidy= ploidy, pheno.file= temp_pheno_file, geno.file=temp_geno_file,
                          format="numeric", n.traits=length(traits), delim=",")
    rm(geno_mat)
    rm(gpoly_df)
    rm(vcf)

  } else {
    stop("Wrong file format")
  }

  data.loco <- set.K(data,LOCO=F,n.core= as.numeric(cores))

  #Delete temp pheno file
  unlink(temp_pheno_file)

  ####Pheno, kinship, PCs from results of GWASpoly
  GE<- data@pheno
  names(GE)
  colnames(GE)[1]<-"Genotype"

  ## kinship
  Kin<- data.loco@K$all

  ## PCs
  PC_all<- eigen(data.loco@K$all)$vectors
  rownames(PC_all) <- rownames(data.loco@K$all)
  PCs<-PC_all[,1:10]
  colnames(PCs)<-c(paste("PC",1:10,sep=""))

  ##
  taxa<-intersect(GE$Genotype,intersect(rownames(PCs),rownames(Kin)))
  PCs<-PCs[which(rownames(PCs) %in% taxa),]
  PCs<-PCs[order(rownames(PCs)),]

  GE<-GE[which(GE$Genotype %in% taxa),]
  GE<-GE[order(GE$Genotype),]

  Kin<-Kin[which(rownames(Kin) %in% taxa),which(rownames(Kin) %in% taxa)] # need check the matrix after this step
  Kin<-Kin[order(rownames(Kin)),order(colnames(Kin))]

  #### calculate BIC
  PC<-as.matrix(PCs)
  K=as.matrix(Kin)

  kin.adj<-posdefmat(K)
  kin.test<-as.matrix(kin.adj)

  for (i in 2:ncol(GE)){

    #model selection
    y=as.numeric(GE[,i])

    BICs<-CalcBIC(y=y,PC=PC,K=kin.test)

    plotBICs<-cbind(rbind.data.frame(BICs$BIC$withK,BICs$BIC$withoutK),rep(c("w/Kinship","no Kinship"),each=nrow(BICs$BIC$withK)))
    colnames(plotBICs)[ncol(plotBICs)]<-"RelationshipMatrix"
    plotBICs$n.PC<-factor(plotBICs$n.PC,levels=c("0","1","2","3","4","5",
                                                 "6","7","8","9","10"))
    plotBICs_kinship <- subset(plotBICs,plotBICs$RelationshipMatrix =="w/Kinship")
    plotBICs_kinship

    p1<-ggplot(plotBICs_kinship, aes(x=n.PC, y=BIC,group=RelationshipMatrix)) +
      geom_line(color="grey")+
      geom_point(shape=21, color="black", fill="#d95f0e", size=3)+
      theme(text=element_text(size=15),axis.text.x = element_text(angle =0),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      labs(x = "Number of PCs",y="BIC")

    ##GWAS based on model selection
    N <- nrow(data@pheno) #Population size
    #Select models depending on ploidy
    if (ploidy > 2) {
      model <- c("additive","1-dom","2-dom","general","diplo-general","diplo-additive")
    }else{
      model <- c("additive", "1-dom")
    }

    BIC_min <- plotBICs_kinship[which.min(plotBICs_kinship$BIC),]
    if(BIC_min$n.PC == 0){params <- set.params(geno.freq = 1 - 5/N)}else{params <- set.params(geno.freq = 1 - 5/N,n.PC = as.numeric(levels(BIC_min$n.PC))[BIC_min$n.PC])}
    data.loco.scan <- GWASpoly(data=data.loco,models=model,traits=colnames(data@pheno[i]),params=params,n.core=as.numeric(cores))
    #Consider adding options for different thresholds
    data2 <- set.threshold(data.loco.scan,method=input$gwas_threshold,level=0.05)

    #Save manhattan plots to list (only for single trait analysis)
    #if length(traits) == 1
    manhattan_plot_list <- list()

    #plot for six models per trait
    manhattan_plot_list[["all"]] <- manhattan.plot(data2,traits=colnames(data@pheno[i]), models = model)+geom_point(size=3)+theme(text = element_text(size = 25),strip.text = element_text(face = "bold"))

    #Output the manhattan plots
    manhattan_plot_list[[input$model_select]]

    #get most significant SNPs per QTL file
    qtl <- get.QTL(data=data2,traits=colnames(data@pheno[i]),bp.window=5e6)
    qtl_d <- data.frame(qtl)

    #get qqplot
    data_qq <- cbind.data.frame(SNP=data.loco.scan@map$Marker,Chr=data.loco.scan@map$Chrom, Pos=data.loco.scan@map$Position,10^(-data.loco.scan@scores[[colnames(data@pheno[i])]]))

    #Save qq_plot info

    CMplot_shiny(data_qq,plot.type="q",col=c(1:8),
                 ylab.pos=2,
                 file.name=colnames(data@pheno[i]),
                 conf.int=FALSE,
                 box=F,multraits=TRUE,file.output=FALSE)

    #plot for each model per trait
    for (j in 1:length(model)) {
      data.loco.scan_2 <- GWASpoly(data=data.loco,models=model[j],
                                   traits=colnames(data@pheno[i]),params=params,n.core= as.numeric(cores))

      data3 <- set.threshold(data.loco.scan_2,method="M.eff",level=0.05)
      manhattan_plot_list[[model[j]]] <- manhattan.plot(data3,traits=colnames(data@pheno[i]))+geom_point(size=3)+theme(text = element_text(size = 25),strip.text = element_text(face = "bold"))
    }

    #Save manhattan plots
    manhattan_plot_list
  }

})
