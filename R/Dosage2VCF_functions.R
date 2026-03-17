#' Helper Functions for Raw MADC Preprocessing in dosage2vcf Module
#'
#' These functions encapsulate the logic for detecting, validating, and preprocessing
#' raw (unprocessed) MADC files before conversion to VCF format.

#' Detect if MADC file is raw (unprocessed)
#'
#' @param madc_path Path to MADC file
#' @return Logical indicating if file is raw (first 7 rows are "*" or blank)
#' @noRd
is_raw_madc <- function(madc_path) {
  lines <- readLines(madc_path, n = 7, warn = FALSE)
  if (length(lines) < 7) return(FALSE)

  first_field <- vapply(strsplit(lines, ",", fixed = TRUE), `[`, character(1), 1)
  first_field <- trimws(first_field)
  all(first_field %in% c("", "*"))
}

#' Validate marker file format
#'
#' @param marker_path Path to marker file
#' @return List with $valid (logical) and $message (character) or $data (validated df)
#' @noRd
read_marker_file <- function(marker_path) {
  marker_df <- read.csv(
    marker_path,
    header = FALSE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(marker_df) == 0) {
    stop("Marker file is empty")
  }

  if (ncol(marker_df) %in% c(3, 4)) {
    header_row <- tolower(trimws(as.character(marker_df[1, ])))
    has_header <- header_row[1] %in% c("cloneid", "clone_id", "clone id") &&
      header_row[2] %in% c("chr", "chromosome", "chrom") &&
      header_row[3] %in% c("pos", "position")

    if (has_header) {
      marker_df <- marker_df[-1, , drop = FALSE]
    }
  }

  if (nrow(marker_df) == 0) {
    stop("Marker file contains only a header row")
  }

  marker_df
}

#' Validate marker file format
#'
#' @param marker_path Path to marker file
#' @return List with $valid (logical) and $message (character) or $data (validated df)
#' @noRd
validate_marker_file <- function(marker_path) {
  tryCatch({
    marker_df <- read_marker_file(marker_path)

    # Check column count
    if (ncol(marker_df) < 3 || ncol(marker_df) > 4) {
      return(list(valid = FALSE,
                  message = "Marker file must have 3 or 4 columns: CloneID, Chr, Pos, [BottomStrand]"))
    }

    # Trim whitespace
    marker_df[,1:3] <- lapply(marker_df[,1:3, drop = FALSE], function(x) trimws(as.character(x)))
    if (ncol(marker_df) == 4) {
      marker_df[,4] <- trimws(as.character(marker_df[,4]))
    }

    if (any(marker_df[,1] == "" | marker_df[,2] == "" | marker_df[,3] == "")) {
      return(list(valid = FALSE,
                  message = "Columns 1-3 cannot contain blank values"))
    }

    # Validate position column (numeric only)
    if (!all(grepl("^[0-9]+$", marker_df[,3]))) {
      return(list(valid = FALSE,
                  message = "Column 3 (Pos) must contain only numeric values"))
    }

    # Validate Chr column (no special characters)
    if (any(grepl("[*#_!.\\-]", marker_df[,2]))) {
      return(list(valid = FALSE,
                  message = "Column 2 (Chr) cannot contain special characters (*#_-!.)"))
    }

    # Check for duplicate CloneIDs
    if (length(unique(marker_df[,1])) != nrow(marker_df)) {
      return(list(valid = FALSE,
                  message = "Column 1 (CloneID) contains duplicate entries"))
    }

    return(list(valid = TRUE, data = marker_df))

  }, error = function(e) {
    return(list(valid = FALSE, message = paste("Error reading file:", e$message)))
  })
}

#' Create CloneID lookup table and generate botloci from 4th column
#'
#' @param marker_df Validated marker data frame
#' @return List with $lookup (named vector), $botloci_ids, and $has_bottom_strand_col
#' @noRd
process_marker_file <- function(marker_df) {
  # Create lookup: old CloneID -> new Chr_Pos format (matching fixMADC)
  lookup <- setNames(
    paste0(marker_df[,2], "_", sprintf("%09d", as.integer(marker_df[,3]))),
    marker_df[,1]
  )

  # Process 4th column if present
  botloci_ids <- NULL
  has_bottom_strand_col <- ncol(marker_df) == 4
  if (has_bottom_strand_col) {
    bottom_strand_mask <- tolower(trimws(as.character(marker_df[,4]))) %in%
      c("y", "yes", "true", "1")

    # Map bottom-strand CloneIDs to the Chr_Pos format expected by BIGr.
    botloci_ids <- unname(lookup[marker_df[bottom_strand_mask, 1]])
  }

  list(
    lookup = lookup,
    botloci_ids = botloci_ids,
    has_bottom_strand_col = has_bottom_strand_col
  )
}

#' Check for markers missing Ref/Alt allele IDs
#'
#' @param madc_path Path to raw MADC file
#' @param marker_cloneids Character vector of CloneIDs from marker file
#' @return Character vector of CloneIDs missing Ref or Alt, or NULL if all present
#' @noRd
check_missing_ref_alt <- function(madc_path, marker_cloneids) {
  raw_madc <- read.csv(madc_path, skip = 7, check.names = FALSE)

  markers_with_ref <- unique(raw_madc$CloneID[grepl("\\|Ref", raw_madc$AlleleID)])
  markers_with_alt <- unique(raw_madc$CloneID[grepl("\\|Alt", raw_madc$AlleleID)])
  markers_complete <- intersect(markers_with_ref, markers_with_alt)

  missing <- setdiff(marker_cloneids, markers_complete)

  if (length(missing) > 0) missing else NULL
}

#' Write normalized marker CSV for BIGr::fixMADC
#'
#' @param marker_df Validated marker data frame
#' @return Path to temporary marker CSV with a header row
#' @noRd
write_marker_file <- function(marker_df) {
  temp_marker <- tempfile(fileext = ".csv")
  marker_to_write <- as.data.frame(marker_df, stringsAsFactors = FALSE)

  colnames(marker_to_write) <- c("CloneID", "Chr", "Pos", "BottomStrand")[seq_len(ncol(marker_to_write))]

  write.csv(marker_to_write, temp_marker, row.names = FALSE, quote = FALSE)
  temp_marker
}

#' Write temporary botloci file
#'
#' @param botloci_ids Character vector of loci on bottom strand
#' @return Path to temporary .botloci file
#' @noRd
write_botloci_file <- function(botloci_ids) {
  temp_botloci <- tempfile(fileext = ".botloci")
  file.create(temp_botloci)

  if (length(botloci_ids) > 0) {
    writeLines(botloci_ids, temp_botloci)
  }

  temp_botloci
}

#' Preprocess raw MADC file with fixMADC
#'
#' @param madc_path Path to raw MADC file
#' @param marker_path Path to marker file
#' @return Path to temporary fixed MADC file
#' @noRd
preprocess_raw_madc <- function(madc_path, marker_path) {
  temp_fixed <- tempfile(fileext = ".csv")

  # Call BIGr::fixMADC
  fixed_df <- BIGr::fixMADC(
    madc.file = madc_path,
    marker.file = marker_path,
    n.summary.columns = NULL,
    output.file = NULL
  )

  # Write to temp file
  write.csv(fixed_df, temp_fixed, row.names = FALSE)

  temp_fixed
}
