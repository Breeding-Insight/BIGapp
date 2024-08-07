context("PCA")

test_that("",{

  # inputs
  input <- list()
  input$dosage_file$datapath <- system.file("iris_DArT_VCF.vcf", package = "BIGapp")
  input$passport_file$datapath <- system.file("iris_passport_file.csv", package = "BIGapp")
  input$pca_ploidy <- 2

  input$group_info <- "Species"
  input$use_cat <- TRUE
  input$color_choice <- "Paired"
  input$pc_X <- "PC1"
  input$pc_Y <- "PC2"
  input$grey_choice <- "Grey"

  all_plots <- pca_data <- data <- list()

  #PCA dropdown
  # Update dropdown menu choices based on uploaded passport file
  info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE)
  info_df[,1] <- as.character(info_df[,1]) #Makes sure that the sample names are characters instead of numeric
  data$info_df <- info_df  # Store info_df in reactiveValues

  # Get selected column name
  selected_col <- input$group_info

  # Extract unique values from the selected column
  unique_values <- unique(data$info_df[[selected_col]])

  #PCA events
  # Get inputs
  geno <- input$dosage_file$datapath
  g_info <- as.character(input$group_info)
  output_name <- input$output_name
  ploidy <- input$pca_ploidy
  PCX <- input$pc_X
  PCY <- input$pc_Y

  #Import genotype information if in VCF format
  vcf <- read.vcfR(geno)

  #Get items in FORMAT column
  info <- vcf@gt[1,"FORMAT"] #Getting the first row FORMAT

  # Apply the function to the first INFO string
  info_ids <- extract_info_ids(info[1])

  #Get the genotype values if the updog dosage calls are present
  if ("UD" %in% info_ids) {
    genomat <- extract.gt(vcf, element = "UD")
    class(genomat) <- "numeric"
  }else{
    #Extract GT and convert to numeric calls
    genomat <- extract.gt(vcf, element = "GT")
    genomat <- apply(genomat, 2, convert_to_dosage)
  }
  rm(vcf) #Remove VCF

  #Start analysis

  # Passport info
  info_df <- read.csv(input$passport_file$datapath, header = TRUE, check.names = FALSE)
  duplicated_samples <- info_df[duplicated(info_df[, 1]), 1]

  # Print the modified dataframe
  row.names(info_df) <- info_df[,1]

  #Plotting
  #First build a relationship matrix using the genotype values
  G.mat.updog <- AGHmatrix::Gmatrix(t(genomat), method = "VanRaden", ploidy = as.numeric(ploidy), missingValue = "NA")

  #PCA
  prin_comp <- prcomp(G.mat.updog, scale = TRUE)
  eig <- factoextra::get_eigenvalue(prin_comp)
  round(sum(eig$variance.percent[1:3]),1)

  ###Simple plots
  # Extract the PC scores
  pc_scores <- prin_comp$x

  # Create a data frame with PC scores
  pc_df <- data.frame(PC1 = pc_scores[, 1], PC2 = pc_scores[, 2],
                      PC3 = pc_scores[, 3], PC4 = pc_scores[, 4],
                      PC5 = pc_scores[, 5], PC6 = pc_scores[, 6],
                      PC7 = pc_scores[, 7], PC8 = pc_scores[, 8],
                      PC9 = pc_scores[, 9], PC10 = pc_scores[, 10])


  # Compute the percentage of variance explained for each PC
  variance_explained <- round(100 * prin_comp$sdev^2 / sum(prin_comp$sdev^2), 1)


  # Retain only samples in common
  row.names(info_df) <- info_df[,1]
  info_df <- info_df[row.names(pc_df),]

  #Add the information for each sample
  pc_df_pop <- merge(pc_df, info_df, by.x = "row.names", by.y = "row.names", all.x = TRUE)


  # Ignore color input if none is entered by user
  pc_df_pop[[g_info]] <- as.factor(pc_df_pop[[g_info]])

  #Update global variable
  pca_dataframes <- pc_df_pop

  # Generate a distinct color palette if g_info is provided
  unique_countries <- unique(pc_df_pop[[g_info]])
  palette <- RColorBrewer::brewer.pal(length(unique_countries), input$color_choice)
  my_palette <- colorRampPalette(palette)(length(unique_countries))

  # Store processed data in reactive values
  pca_data$pc_df_pop <- pc_df_pop
  pca_data$variance_explained <- variance_explained
  pca_data$my_palette <- my_palette

  ##2D PCA plotting
  # Generate colors
  unique_countries <- unique(pca_data$pc_df_pop[[input$group_info]])
  palette <- RColorBrewer::brewer.pal(length(unique_countries), input$color_choice)
  my_palette <- colorRampPalette(palette)(length(unique_countries))



  # Define a named vector to map input labels to grey values
  label_to_value <- c("Light Grey" = "grey80",
                      "Grey" = "grey60",
                      "Dark Grey" = "grey40",
                      "Black" = "black")

  # Get the corresponding value based on the selected grey
  selected_grey <- label_to_value[[input$grey_choice]]

  #Set factor
  if (!input$use_cat && is.null(my_palette)) {
    print("No Color Info")
  }else{
    pca_data$pc_df_pop[[input$group_info]] <- as.factor(pca_data$pc_df_pop[[input$group_info]])
  }

  input$cat_color <- "versicolor"
  # Similar plotting logic here

  cat_colors <- c(input$cat_color, "grey")
  plot <- {if(!is.null(input$group_info) & input$group_info != "")
    ggplot(pca_data$pc_df_pop, aes(x = pca_data$pc_df_pop[[input$pc_X]],
                                   y = pca_data$pc_df_pop[[input$pc_Y]],
                                   color = factor(pca_data$pc_df_pop[[input$group_info]]))) else
                                     ggplot(pca_data$pc_df_pop, aes(x = pca_data$pc_df_pop[[input$pc_X]],
                                                                    y = pca_data$pc_df_pop[[input$pc_Y]]))} +
    geom_point(size = 2, alpha = 0.8) +
    {if(input$use_cat) scale_color_manual(values = setNames(c(my_palette, "grey"), cat_colors), na.value = selected_grey) else
      if(!is.null(my_palette)) scale_color_manual(values = my_palette)} +
    guides(color = guide_legend(override.aes = list(size = 5.5), nrow = 17)) +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 16)
    ) +
    labs(
      x = paste0(input$pc_X, "(", pca_data$variance_explained[as.numeric(substr(input$pc_X, 3, 3))], "%)"),
      y = paste0(input$pc_Y, "(", pca_data$variance_explained[as.numeric(substr(input$pc_Y, 3, 3))], "%)"),
      color = input$group_info
    )

  plot  # Assign the plot to your reactiveValues

  #3D PCA plotting
  #Generate colors
  unique_countries <- unique(pca_data$pc_df_pop[[input$group_info]])
  palette <- brewer.pal(length(unique_countries),input$color_choice)
  my_palette <- colorRampPalette(palette)(length(unique_countries))

  tit = paste0('Total Explained Variance =', sum(pca_data$variance_explained[1:3]))

  fig <- plot_ly(pca_data$pc_df_pop, x = ~PC1, y = ~PC2, z = ~PC3, color = pca_data$pc_df_pop[[input$group_info]],
                 colors = my_palette) %>%
    add_markers(size = 12, text = paste0("Sample:",pca_data$pc_df_pop$Row.names))

  fig <- fig %>%
    layout(
      title = tit,
      scene = list(bgcolor = "white")
    )

  fig # Return the Plotly object here

  #PCA scree plot
  var_explained <- pca_data$variance_explained

  # Create a data frame for plotting
  plot_data <- data.frame(PC = 1:10, Variance_Explained = var_explained[1:10])

  # Use ggplot for plotting
  plot <- ggplot(plot_data, aes(x = PC, y = Variance_Explained)) +
    geom_bar(stat = "identity", fill = "lightblue", alpha = 0.9, color = "black") +  # Bars with some transparency
    geom_line(color = "black") +  # Connect points with a line
    geom_point(color = "black") +  # Add points on top of the line for emphasis
    scale_x_continuous(breaks = 1:10, limits = c(0.5, 10.5)) +
    xlab("Principal Component") +
    ylab("% Variance Explained") +
    ylim(0, 100) +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 16)
    )

  plot
})
