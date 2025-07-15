#' Perform a Sanity Check on a VCF File
#'
#' This function performs a series of checks on a VCF file to ensure its validity and integrity. It verifies the presence of required headers, columns, and data fields, and checks for common issues such as missing or malformed data.
#'
#' @param vcf_path A character string specifying the path to the VCF file. The file can be plain text or gzipped.
#' @param n_data_lines An integer specifying the number of data lines to sample for detailed checks. Default is 100.
#' @param depth_support_fields (Optional) A character vector of fields that are expected to be present in the FORMAT column for allele counts. Default is `c("AD", "RA", "AO", "RO", "NR", "NV", "SB", "F1R2", "F2R1")`.
#' @param max_markers An integer specifying the maximum number of markers allowed in the VCF file. Default is 10,000.
#' @param verbose A logical value indicating whether to print detailed messages during the checks. Default is FALSE.
#'
#' @return A list containing:
#' - `checks`: A named vector indicating the results of each check (TRUE or FALSE).
#' - `messages`: A data frame containing messages for each check, indicating success or failure.
#' - `duplicates`: A list containing any duplicated sample or marker IDs found in the VCF file.
#' - `ploidy_max`: The maximum ploidy detected from the genotype field, if applicable.
#'
#' @details The function performs the following checks:
#' - **VCF_header**: Verifies the presence of the `##fileformat` header.
#' - **VCF_compressed**: Checks if the VCF file is compressed and if the extension is correct.
#' - **VCF_columns**: Ensures required columns (`#CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`) are present.
#' - **max_markers**: Checks if the total number of markers exceeds the specified limit.
#' - **unique_FORMAT**: Ensures that the FORMAT fields are consistent across sampled markers.
#' - **GT**: Verifies the presence of the `GT` (genotype) field in the FORMAT column.
#' - **allele_counts**: Checks for allele-level count fields (e.g., `AD`, `RA`, `AO`, `RO`).
#' - **samples**: Ensures sample/genotype columns are present.
#' - **chrom_info** and **pos_info**: Verifies the presence of `CHROM` and `POS` columns.
#' - **ref_alt**: Ensures `REF` and `ALT` fields contain valid nucleotide codes.
#' - **multiallelics**: Identifies multiallelic sites (ALT field with commas).
#' - **phased_GT**: Checks for phased genotypes (presence of `|` in the `GT` field).
#' - **duplicated_samples**: Checks for duplicated sample IDs.
#' - **duplicated_markers**: Checks for duplicated marker IDs.
#'
#' @importFrom stats setNames
#'
#' @export
vcf_sanity_check <- function(
    vcf_path,
    n_data_lines = 100,
    max_markers = 10000,
    depth_support_fields = c("AD", "RA", "AO", "RO", "NR", "NV", "SB", "F1R2", "F2R1"),
    verbose = FALSE) {
  
  # --- Prepare result vector ---
  checks_names <- c(
    "VCF_header",
    "VCF_compressed",
    "VCF_columns",
    "max_markers",
    "unique_FORMAT",
    "GT",
    "allele_counts",
    "samples",
    "chrom_info",
    "pos_info",
    "ref_alt",
    "multiallelics",
    "phased_GT",
    "duplicated_samples",
    "duplicated_markers",
    "mixed_ploidies"
  )
  checks <- setNames(rep(NA, length(checks_names)), checks_names)
  
  if (!file.exists(vcf_path)) stop("File does not exist.")
  
  is_gz_ext <- grepl("\\.gz$", vcf_path)
  is_gz <- is_compressed_file(vcf_path)
  
  con <- if (is_gz == "gzip (.gz)") {
    checks["VCF_compressed"] <- TRUE
    if(!is_gz_ext) if(verbose) warning("File is compressed with gzip (.gz), but does not have .gz extension.")
    gzfile(vcf_path, open = "rt") 
  } else if(is_gz == "bzip2 (.bz2)") {
    if (verbose) warning("File is compressed th bzip2 (.bz2), which is not supported.")
    checks["VCF_compressed"] <- FALSE
  } else if (is_gz == "xz (.xz)") {
    if (verbose) warning("File is compressed with xz (.xz), which is not supported.")
    checks["VCF_compressed"] <- FALSE
  } else {
    checks["VCF_compressed"] <- FALSE
    file(vcf_path, open = "r")
  }
  
  lines <- readLines(con, warn = FALSE)
  close(con)
  
  # Container for duplicated IDs
  duplicates <- list(
    duplicated_samples = character(0),
    duplicated_markers = character(0)
  )
  
  # --- Header checks ---
  header_lines <- grep("^##", lines, value = TRUE)
  if (!any(grepl("^##fileformat=VCFv", header_lines))) {
    checks["VCF_header"] <- FALSE
    if (verbose) warning("Missing ##fileformat header.")
  } else {
    checks["VCF_header"] <- TRUE
    if (verbose) cat("VCF header is present.\n")
  }
  
  # --- Column header line ---
  column_header_line <- grep("^#CHROM", lines, value = TRUE)
  if (length(column_header_line) != 1) stop("Missing or multiple #CHROM lines.")
  column_names <- unlist(strsplit(column_header_line, "\t"))
  has_genotypes <- length(column_names) > 8
  
  required_columns <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  checks["VCF_columns"] <- all(required_columns %in% column_names[1:8])
  if (checks["VCF_columns"]) {
    if (verbose) cat("Required VCF columns are present.\n")
  } else {
    if (verbose) warning("Missing one or more required VCF columns.")
  }
  
  # --- Total marker count ---
  data_line_indices <- grep("^[^#]", lines)
  total_markers <- length(data_line_indices)
  if (verbose) cat(sprintf("Total markers (data rows): %d\n", total_markers))
  
  checks["max_markers"] <- total_markers <= max_markers
  if (!checks["max_markers"]) {
    warning(sprintf("More than %d markers found. Consider subsampling.", max_markers))
  }
  
  # --- Check for duplicated marker IDs ---
  if (total_markers > 0) {
    marker_ids <- sapply(lines[data_line_indices], function(line) {
      fields <- strsplit(line, "\t")[[1]]
      if (length(fields) >= 3) fields[3] else NA
    })
    marker_ids <- marker_ids[!is.na(marker_ids)]
    duplicated_markers <- marker_ids[duplicated(marker_ids)]
    duplicates$duplicated_markers <- duplicated_markers
    if (length(duplicated_markers) > 0) {
      checks["duplicated_markers"] <- TRUE
      if (verbose) warning("Duplicated marker IDs found: ", paste(head(duplicated_markers, 10), collapse = ", "), "...")
    } else {
      if (verbose) cat("No duplicated marker IDs.\n")
      checks["duplicated_markers"] <- FALSE
    }
  }
  
  # --- FORMAT field checks (GT, AD, etc.) ---
  if (has_genotypes) {
    checks["samples"] <- TRUE
    sample_indices <- sample(data_line_indices, n_data_lines)
    format_fields <- character()
    
    chrom_pos <- list() # this list will store the CHROM and POS for n_data_lines sample markers
    format_keys <- list() # this list store the format field for n_data_lines sample markers
    for (i in seq_along(sample_indices)) {
      fields <- strsplit(lines[sample_indices[i]], "\t")[[1]]
      if (length(fields) >= 9) {
        format_keys[[i]] <- unlist(strsplit(fields[9], ":"))
        chrom_pos[[i]] <- fields[1:2]
      }
    }
    
    # Check if all markers sampled format are identical with the first
    checks["unique_FORMAT"] <- all(sapply(format_keys[-1], function(x) identical(x, format_keys[[1]])))
    
    # --- CHROM and POS column checks ---
    chrom_pos <- do.call(rbind, chrom_pos)
    checks["chrom_info"] <- all(chrom_pos[,1] != "." | chrom_pos[,1] != "" | !is.na(chrom_pos[,1]))
    checks["pos_info"] <- all(chrom_pos[,2] != "." | chrom_pos[,2] != "" | !is.na(chrom_pos[,2]))
    
    if (checks["chrom_info"] && checks["pos_info"]) {
      if (verbose) cat("Both CHROM and POS information are present.\n")
    } else {
      if (!checks["chrom_info"]) warning(" 'CHROM' information is missing for at least one marker.")
      if (!checks["pos_info"]) warning(" 'POS' information is missing for at least one marker.")
    }
    
    format_fields <- unique(unlist(format_keys))
    
    # GT check
    checks["GT"] <- "GT" %in% format_fields
    if (verbose) {
      if (checks["GT"]) cat("FORMAT field 'GT' (genotype) is present.\n") else warning("FORMAT field 'GT' is missing.")
    }
    # Allele counts check
    checks["allele_counts"] <- any(depth_support_fields %in% format_fields)
    if (checks["allele_counts"]) {
      if (verbose) {
        cat(sprintf(
          "Allele count FORMAT field(s) found: %s\n",
          paste(intersect(depth_support_fields, format_fields), collapse = ", ")
        ))
      }
    } else {
      warning(paste(" No required allele-level count fields found. e.g.", depth_support_fields))
    }
    
    # Optional: phased GT (presence of '|' instead of '/' in genotypes)
    phased_lines <- grep("\\|", lines[sample_indices])
    checks["phased_GT"] <- length(phased_lines) > 0
    
    # --- Check for duplicated sample names ---
    sample_names <- column_names[10:length(column_names)]
    duplicated_samples <- sample_names[duplicated(sample_names)]
    duplicates$duplicated_samples <- duplicated_samples
    if (length(duplicated_samples) > 0) {
      if (verbose) warning("Duplicated sample names found: ", paste(duplicated_samples, collapse = ", "))
      checks["duplicated_samples"] <- TRUE
    } else {
      if (verbose) cat("No duplicated sample names.\n")
      checks["duplicated_samples"] <- FALSE
    }
  } else {
    checks["samples"] <- FALSE
    checks["GT"] <- FALSE
    checks["allele_counts"] <- FALSE
    checks["phased_GT"] <- FALSE
    checks["duplicated_samples"] <- FALSE
    
    warning("No sample/genotype columns found.")
  }
  
  # --- Ploidy inference (based on GT) ---
  ploidy_max <- NA # default if GT not found or no valid genotypes
  
  if (checks["GT"]) {
    ploidy_values <- c()
    
    for (i in seq_along(sample_indices)) {
      fields <- strsplit(lines[sample_indices[i]], "\t")[[1]]
      
      # Skip if not enough fields
      if (length(fields) < 10) next
      
      format_keys <- unlist(strsplit(fields[9], ":"))
      gt_index <- which(format_keys == "GT")
      
      # Loop through samples (from column 10 onward)
      for (sample_field in fields[10:length(fields)]) {
        sample_fields <- unlist(strsplit(sample_field, ":"))
        if (length(sample_fields) >= gt_index) {
          gt_raw <- sample_fields[gt_index]
          if (grepl("[/|]", gt_raw)) {
            ploidy <- length(unlist(strsplit(gt_raw, "[/|]")))
            ploidy_values <- c(ploidy_values, ploidy)
          }
        }
      }
    }
    
    if (length(ploidy_values) > 0) {
      ploidy_max <- max(ploidy_values, na.rm = TRUE)
      if (verbose) cat("Highest ploidy detected from GT field:", ploidy_max, "\n")
    }
    if (length(unique(ploidy_values)) > 1) {
      checks["mixed_ploidies"] <- TRUE
      if (verbose) cat("Mixed ploidies detected\n")
    } else {
      checks["mixed_ploidies"] <- FALSE
    }
  }

  # --- REF/ALT basic check on sample rows ---
  sample_lines <- lines[head(data_line_indices, n_data_lines)]
  ref_alt_valid <- sapply(sample_lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) >= 5) {
      ref <- fields[4]
      alt <- fields[5]
      grepl("^[ACGTN]+$", ref) && grepl("^[ACGTN.,<>]+$", alt)
    } else {
      FALSE
    }
  })
  checks["ref_alt"] <- all(ref_alt_valid)
  
  # --- Multiallelic site check (ALT with ',' separator) ---
  multiallelic_flags <- grepl(",", sapply(sample_lines, function(line) strsplit(line, "\t")[[1]][5]))
  checks["multiallelics"] <- any(multiallelic_flags)
  
  # --- Compile messages ---
  
  # Messages in case of failure
  messages <- data.frame(
    "VCF_header" = c(
      "VCF header is missing. Please check the file format\n",
      "VCF header is present\n"
    ),
    "VCF_compressed" = c(
      "VCF is uncompressed or compressed with non-supported format\n",
      "VCF is .gz compressed. Check if has the proper extension\n"
    ),
    "VCF_columns" = c(
      "Required VCF columns are missing. Please check the file format\n",
      "Required VCF columns are present\n"
    ),
    "max_markers" = c(
      "More than 10,000 markers found. Consider subsampling or running in HPC\n",
      "Less than maximum number of markers found\n"
    ),
    "unique_FORMAT" = c(
      "FORMAT fields are not consistent across sampled markers\n",
      "FORMAT fields are consistent across sampled markers\n"
    ),
    "GT" = c(
      "Genotype information is not available in the VCF file\n",
      "Genotype information is available in the VCF file\n"
    ),
    "allele_counts" = c(
      "Required field for allele counts are not available in the VCF file\n",
      "Required field for allele counts are available in the VCF file\n"
    ),
    "samples" = c(
      "Sample information is not available in the VCF file\n",
      "Sample information is available in the VCF file\n"
    ),
    "chrom_info" = c(
      "Chromosome information is not available in the VCF file\n",
      "Chromosome information is available in the VCF file\n"
    ),
    "pos_info" = c(
      "Position information is not available in the VCF file\n",
      "Position information is available in the VCF file\n"
    ),
    "ref_alt" = c(
      "REF/ALT fields contain invalid nucleotide codes\n",
      "REF/ALT fields are valid\n"
    ),
    "multiallelics" = c(
      "Multiallelic sites not found in the VCF file\n",
      "Multiallelic sites found in the VCF file\n"
    ),
    "phased_GT" = c(
      "Phased genotypes (|) are not present in the VCF file\n",
      "Phased genotypes (|) are present in the VCF file\n"
    ),
    "duplicated_samples" = c(
      "No duplicated sample IDs found",
      paste("Duplicated sample IDs found: ", paste(duplicates$duplicated_samples, collapse = ", "))
    ),
    "duplicated_markers" = c(
      "No duplicated marker IDs found",
      paste("Duplicated marker IDs found: ", paste(duplicates$duplicated_markers, collapse = ", "))
    ),
    "mixed_ploidies" = c(
      "No mixed ploidies detected",
      "Mixed ploidies detected"
    )
  )
  rownames(messages) <- c("false", "true")
  
  # --- Done ---
  if (verbose) cat("Sanity check complete.\n")
  return(structure(
    list(
      checks = checks, messages = messages,
      duplicates = duplicates, ploidy_max = ploidy_max
    ),
    class = "vcf_sanity_check"
  ))
}

#' Generate Sanity Check Messages for VCF Files
#'
#' This function generates messages based on the results of a VCF sanity check.
#'
#' @param checks A `vcf_sanity_check` object containing the results of the sanity check.
#' @param error_if_false (Optional) A character vector of checks that must be TRUE for the VCF file to pass validation.
#' @param error_if_true (Optional) A character vector of checks that must be FALSE for the VCF file to pass validation.
#' @param warning_if_true (Optional) A character vector of checks that should trigger a warning if TRUE.
#' @param warning_if_false (Optional) A character vector of checks that should trigger a warning if FALSE.
#' @param input_ploidy (Optional) An integer specifying the expected ploidy of the VCF file. If provided, the function will compare it with the detected ploidy.
#'
#' @details This function uses the results of a VCF sanity check to display appropriate messages to the user. It checks for required conditions that must be met (TRUE or FALSE) and optionally validates the ploidy of the VCF file. If any issues are detected, the function displays alerts using `shinyalert`.
#'
#' @return A logical value indicating whether any of the specified checks failed (TRUE) or passed (FALSE). If any required checks are not met, an error alert is displayed. If any warnings are triggered, a warning alert is displayed.
#'
#' @importFrom shinyalert shinyalert
#'
#' @export
vcf_sanity_messages <- function(
    checks,
    error_if_false = c(
      "VCF_header",
      "VCF_columns",
      "max_markers",
      "GT",
      "allele_counts",
      "samples",
      "chrom_info",
      "pos_info",
      "ref_alt",
      "mixed_ploidies"
    ), 
    error_if_true = c( "phased_GT",
                        "duplicated_samples",
                        "duplicated_markers"),
    warning_if_true = c("multiallelics"),
    warning_if_false = c("unique_FORMAT"),
    input_ploidy = NULL) {
  
  # check must be from vcf_sanity_check
  if (!inherits(checks, "vcf_sanity_check")) {
    stop("check must be a vcf_sanity_check object.")
  }
  
  it_should_brake <- vector()
  # Check required TRUE
  if (length(error_if_false) > 0 && !all(checks$checks[error_if_false])) {
    it_should_brake <- c(it_should_brake, !checks$checks[error_if_false])
    shinyalert(
      title = "File Error",
      text = paste(checks$message[1, error_if_false][which(!checks$checks[error_if_false])],
                   collapse = ";"
      ),
      size = "xs",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = TRUE,
      type = "error",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#004192",
      showCancelButton = FALSE,
      imageUrl = "",
      animation = TRUE,
    )
  }
  
  # Check warning TRUE
  if (length(warning_if_true) > 0 && all(checks$checks[warning_if_true])) {
    shinyalert(
      title = "Warning",
      text = paste(checks$message[2, warning_if_true][which(checks$checks[warning_if_true])],
                   collapse = ";"
      ),
      size = "xs",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = TRUE,
      type = "warning",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#004192",
      showCancelButton = FALSE,
      imageUrl = "",
      animation = TRUE,
    )
  }
  
  # Check required FALSE
  if (length(error_if_true) > 0 && any(checks$checks[error_if_true])) {
    it_should_brake <- c(it_should_brake, checks$checks[error_if_true])
    
    shinyalert(
      title = "File Error",
      text = paste(checks$message[2, error_if_true][which(checks$checks[error_if_true])],
                   collapse = ";"
      ),
      size = "xs",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = TRUE,
      type = "error",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#004192",
      showCancelButton = FALSE,
      imageUrl = "",
      animation = TRUE,
    )
  }
  
  # Check warning FALSE
  if (length(warning_if_false) > 0 && !all(checks$checks[warning_if_false])) {
    shinyalert(
      title = "Warning",
      text = paste(checks$message[1, warning_if_false][which(!checks$checks[warning_if_false])],
                   collapse = ";"
      ),
      size = "xs",
      closeOnEsc = TRUE,
      closeOnClickOutside = FALSE,
      html = TRUE,
      type = "warning",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      confirmButtonCol = "#004192",
      showCancelButton = FALSE,
      imageUrl = "",
      animation = TRUE,
    )
  }
  
  if (!is.null(input_ploidy)) {
    it_should_brake <- c(it_should_brake, checks$ploidy_max != input_ploidy)
    
    if (checks$ploidy_max != input_ploidy) {
      shinyalert(
        title = "File Error",
        text = "Informed ploidy doesn't match genotypes in the VCF file",
        size = "xs",
        closeOnEsc = TRUE,
        closeOnClickOutside = FALSE,
        html = TRUE,
        type = "info",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
        confirmButtonCol = "#004192",
        showCancelButton = FALSE,
        imageUrl = "",
        animation = TRUE,
      )
    }
  }
  return(any(it_should_brake))
}

#' Print Method for vcf_sanity_check Objects
#'
#' Displays a summary table of the VCF sanity check results, showing the check name, result (TRUE/FALSE), and the corresponding message for each selected check.
#'
#' @param x A vcf_sanity_check object as returned by `vcf_sanity_check()`.
#' @param checks A character vector of check names to display. Defaults to all checks in the object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return The input object, invisibly. Called for its side effect (printing to console).
#' 
#' @method print vcf_sanity_check
#'
#' @export
print.vcf_sanity_check <- function(x, checks = names(x$checks), ...) {
  if (!inherits(x, "vcf_sanity_check")) stop("Object must be of class 'vcf_sanity_check'.")
  
  # Prepare a data.frame with check name, result, and message
  res <- data.frame(
    Check = checks,
    Result = as.logical(x$checks[checks]),
    Message = mapply(
      function(chk, val) {
        # TRUE is row 'true', FALSE is row 'false'
        row_idx <- ifelse(isTRUE(val), "true", "false")
        as.character(x$messages[row_idx, chk])
      },
      chk = checks,
      val = x$checks[checks],
      SIMPLIFY = TRUE
    ),
    stringsAsFactors = FALSE
  )
  print(res, row.names = FALSE)
  invisible(x)
}
