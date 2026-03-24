# Load required libraries
library(Biostrings)
library(dplyr)
library(readr)

# Function to generate realistic count data with many zeros
generate_realistic_counts <- function(n_samples, allele_type = "main", min_count = 0, max_count = 200) {
  counts <- numeric(n_samples)

  if (allele_type == "main") {
    # Main alleles (Ref_0001, Alt_0002): 70% chance of being observed in a sample
    for (i in 1:n_samples) {
      if (runif(1) < 0.7) {  # 70% chance of non-zero
        counts[i] <- sample(1:max_count, 1)
      } else {
        counts[i] <- 0
      }
    }
  } else {
    # Match alleles (RefMatch, AltMatch): 30% chance of being observed in a sample
    for (i in 1:n_samples) {
      if (runif(1) < 0.3) {  # 30% chance of non-zero
        counts[i] <- sample(1:100, 1)  # Lower max count for matches
      } else {
        counts[i] <- 0
      }
    }
  }

  return(counts)
}

# Function to generate simulated allele data using FASTA and SNP info
generate_allele_data <- function(fasta_file, snp_file, n_clones = 1000, n_samples = 100, force_include = NULL) {

  # Read FASTA file
  cat("Reading FASTA file...\n")
  sequences <- readDNAStringSet(fasta_file)
  seq_names <- names(sequences)
  seq_chars <- as.character(sequences)
  names(seq_chars) <- seq_names

  # Read SNP information
  cat("Reading SNP information file...\n")
  snp_data <- read_csv(snp_file)
  cat("Found", nrow(snp_data), "SNPs in the file\n")

  # Handle forced inclusions
  selected_snps <- data.frame()

  if (!is.null(force_include) && length(force_include) > 0) {
    cat("Forcing inclusion of", length(force_include), "specific CloneIDs...\n")

    # Try to match force_include with either Panel_markerID or BI_markerID
    forced_snps <- snp_data[snp_data$Panel_markerID %in% force_include |
                              snp_data$BI_markerID %in% force_include, ]

    if (nrow(forced_snps) > 0) {
      cat("Found", nrow(forced_snps), "of", length(force_include), "requested CloneIDs in the data\n")
      selected_snps <- forced_snps

      # Show which ones were found
      found_ids <- unique(c(forced_snps$Panel_markerID, forced_snps$BI_markerID))
      found_requested <- intersect(force_include, found_ids)
      missing_requested <- setdiff(force_include, found_ids)

      if (length(missing_requested) > 0) {
        cat("Warning: Could not find these requested CloneIDs:", paste(missing_requested, collapse = ", "), "\n")
      }
    } else {
      cat("Warning: None of the requested CloneIDs were found in the data\n")
    }
  }

  # Calculate how many more SNPs we need to sample randomly
  remaining_needed <- n_clones - nrow(selected_snps)

  if (remaining_needed > 0) {
    # Get remaining SNPs (exclude already selected ones)
    if (nrow(selected_snps) > 0) {
      remaining_snps <- snp_data[!(snp_data$Panel_markerID %in% selected_snps$Panel_markerID |
                                     snp_data$BI_markerID %in% selected_snps$BI_markerID), ]
    } else {
      remaining_snps <- snp_data
    }

    # Randomly sample from remaining SNPs
    if (remaining_needed > nrow(remaining_snps)) {
      cat("Warning: Requested", remaining_needed, "additional CloneIDs but only", nrow(remaining_snps), "remaining SNPs available.\n")
      cat("Using all remaining SNPs.\n")
      additional_snps <- remaining_snps
    } else {
      cat("Randomly selecting", remaining_needed, "additional SNPs from", nrow(remaining_snps), "remaining.\n")
      selected_indices <- sample(1:nrow(remaining_snps), remaining_needed, replace = FALSE)
      additional_snps <- remaining_snps[selected_indices, ]
    }

    # Combine forced and random selections
    selected_snps <- rbind(selected_snps, additional_snps)
  } else if (remaining_needed < 0) {
    cat("Warning: More forced CloneIDs (", nrow(selected_snps), ") than requested total (", n_clones, ")\n")
    cat("Using only the forced CloneIDs\n")
    # Keep only the forced ones, no random sampling
  }

  cat("Final selection: ", nrow(selected_snps), "CloneIDs\n")

  # Generate sample names
  sample_names <- paste0("Sample_", sprintf("%03d", 1:n_samples))

  # Initialize results dataframe
  results <- data.frame()

  cat("Generating allele data for", nrow(selected_snps), "CloneIDs...\n")

  for (i in 1:nrow(selected_snps)) {
    if (i %% 100 == 0) cat("Processing CloneID", i, "of", nrow(selected_snps), "\n")

    # Get SNP information
    panel_id <- selected_snps$Panel_markerID[i]
    bi_id <- selected_snps$BI_markerID[i]
    chr <- selected_snps$Chr[i]
    pos <- selected_snps$Pos[i]
    ref_base <- toupper(selected_snps$Ref[i])
    alt_base <- toupper(selected_snps$Alt[i])

    # Use BI_markerID as CloneID (or Panel_markerID if you prefer)
    clone_id <- bi_id

    # Find Ref and Alt sequences directly from FASTA (use as-is, no modification)
    ref_keys <- c(paste0(bi_id, "|Ref_0001"), paste0(panel_id, "|Ref_0001"))
    alt_keys <- c(paste0(bi_id, "|Alt_0002"), paste0(panel_id, "|Alt_0002"))

    ref_sequence <- NULL
    alt_sequence <- NULL

    # Exact match for Ref
    for (key in ref_keys) {
      if (key %in% seq_names) { ref_sequence <- seq_chars[[key]]; break }
      matches <- grep(key, seq_names, value = TRUE, fixed = TRUE)
      if (length(matches) > 0) { ref_sequence <- seq_chars[[matches[1]]]; break }
    }

    # Exact match for Alt
    for (key in alt_keys) {
      if (key %in% seq_names) { alt_sequence <- seq_chars[[key]]; break }
      matches <- grep(key, seq_names, value = TRUE, fixed = TRUE)
      if (length(matches) > 0) { alt_sequence <- seq_chars[[matches[1]]]; break }
    }

    # Fallback: if Ref/Alt not found separately, derive from a base sequence
    if (is.null(ref_sequence) || is.null(alt_sequence)) {
      possible_names <- c(panel_id, bi_id,
                          paste0(chr, "_", pos),
                          paste0("chr", chr, "_", pos),
                          paste0(chr, ":", pos))
      base_sequence <- NULL
      seq_found <- FALSE
      for (name in possible_names) {
        if (name %in% seq_names) { base_sequence <- seq_chars[[name]]; seq_found <- TRUE; break }
      }
      if (!seq_found) {
        for (name in possible_names) {
          matches <- grep(name, seq_names, value = TRUE)
          if (length(matches) > 0) { base_sequence <- seq_chars[[matches[1]]]; seq_found <- TRUE; break }
        }
      }
      if (!seq_found) {
        if (length(seq_chars) > 0) {
          base_sequence <- seq_chars[[((i-1) %% length(seq_chars)) + 1]]
        } else {
          seq_length <- sample(70:80, 1)
          base_sequence <- paste(sample(c("A", "T", "G", "C"), seq_length, replace = TRUE), collapse = "")
        }
      }
      # Trim and adjust local position for fallback sequences
      if (nchar(base_sequence) > 80) {
        if (pos <= nchar(base_sequence)) {
          start_pos <- max(1, pos - 40)
          end_pos   <- min(nchar(base_sequence), start_pos + 79)
          start_pos <- max(1, end_pos - 79)
          base_sequence <- substr(base_sequence, start_pos, end_pos)
          pos <- pos - start_pos + 1
        } else {
          base_sequence <- substr(base_sequence, 1, 80)
          pos <- sample(20:60, 1)
        }
      }
      if (pos > nchar(base_sequence) || pos < 1) pos <- sample(1:nchar(base_sequence), 1)
      ref_sequence <- base_sequence
      substr(ref_sequence, pos, pos) <- ref_base
      alt_sequence <- ref_sequence
      substr(alt_sequence, pos, pos) <- alt_base
    }

    # Determine local SNP position (used only to guide RefMatch/AltMatch SNP placement)
    # Compare Ref and Alt to find the differing position
    ref_chars_cmp <- strsplit(toupper(ref_sequence), "")[[1]]
    alt_chars_cmp <- strsplit(toupper(alt_sequence), "")[[1]]
    diff_pos <- which(ref_chars_cmp != alt_chars_cmp)
    local_snp_pos <- if (length(diff_pos) > 0) diff_pos[1] else ceiling(nchar(ref_sequence) / 2)

    # Generate realistic observation counts for samples (with many zeros)
    ref_counts <- generate_realistic_counts(n_samples, "main")
    alt_counts <- generate_realistic_counts(n_samples, "main")

    # Create Ref_0001 row
    ref_row <- data.frame(
      AlleleID = paste0(clone_id, "|Ref_0001"),
      CloneID = clone_id,
      AlleleSequence = ref_sequence,
      stringsAsFactors = FALSE
    )
    # Add sample counts
    for (j in 1:n_samples) {
      ref_row[[sample_names[j]]] <- ref_counts[j]
    }

    # Create Alt_0002 row
    alt_row <- data.frame(
      AlleleID = paste0(clone_id, "|Alt_0002"),
      CloneID = clone_id,
      AlleleSequence = alt_sequence,
      stringsAsFactors = FALSE
    )
    # Add sample counts
    for (j in 1:n_samples) {
      alt_row[[sample_names[j]]] <- alt_counts[j]
    }

    # Add to results
    results <- rbind(results, ref_row, alt_row)

    # Generate RefMatch entries (0-3 random)
    n_ref_matches <- sample(0:3, 1)
    if (n_ref_matches > 0) {
      for (k in 1:n_ref_matches) {
        # RefMatch: introduce a random SNP at a non-SNP position (Ref AlleleSequence unchanged)
        available_positions <- setdiff(1:nchar(ref_sequence), local_snp_pos)
        if (length(available_positions) > 0) {
          ref_match_pos <- sample(available_positions, 1)
          ref_match_base <- substr(ref_sequence, ref_match_pos, ref_match_pos)
          new_bases <- setdiff(c("A", "T", "G", "C"), ref_match_base)
          new_base <- sample(new_bases, 1)

          ref_match_seq <- ref_sequence
          substr(ref_match_seq, ref_match_pos, ref_match_pos) <- new_base
        } else {
          ref_match_seq <- ref_sequence
        }

        # RefMatch sequences are rarer - more zeros
        ref_match_counts <- generate_realistic_counts(n_samples, "match")

        ref_match_row <- data.frame(
          AlleleID = paste0(clone_id, "|RefMatch_", sprintf("%04d", k)),
          CloneID = clone_id,
          AlleleSequence = ref_match_seq,
          stringsAsFactors = FALSE
        )
        # Add sample counts
        for (j in 1:n_samples) {
          ref_match_row[[sample_names[j]]] <- ref_match_counts[j]
        }

        results <- rbind(results, ref_match_row)
      }
    }

    # Generate AltMatch entries (0-3 random)
    n_alt_matches <- sample(0:3, 1)
    if (n_alt_matches > 0) {
      for (k in 1:n_alt_matches) {
        # AltMatch: introduce a random SNP at a non-SNP position (Alt AlleleSequence unchanged)
        available_positions <- setdiff(1:nchar(alt_sequence), local_snp_pos)
        if (length(available_positions) > 0) {
          alt_match_pos <- sample(available_positions, 1)
          alt_match_base <- substr(alt_sequence, alt_match_pos, alt_match_pos)
          new_bases <- setdiff(c("A", "T", "G", "C"), alt_match_base)
          new_base <- sample(new_bases, 1)

          alt_match_seq <- alt_sequence
          substr(alt_match_seq, alt_match_pos, alt_match_pos) <- new_base
        } else {
          alt_match_seq <- alt_sequence
        }

        # AltMatch sequences are rarer - more zeros
        alt_match_counts <- generate_realistic_counts(n_samples, "match")

        alt_match_row <- data.frame(
          AlleleID = paste0(clone_id, "|AltMatch_", sprintf("%04d", k)),
          CloneID = clone_id,
          AlleleSequence = alt_match_seq,
          stringsAsFactors = FALSE
        )
        # Add sample counts
        for (j in 1:n_samples) {
          alt_match_row[[sample_names[j]]] <- alt_match_counts[j]
        }

        results <- rbind(results, alt_match_row)
      }
    }
  }

  return(results)
}
# Simple function to add random lowercase bases to sequences
add_random_lowercase <- function(data, lowercase_prob = 0.15) {
  # Function to convert random bases to lowercase in a single sequence
  convert_sequence <- function(sequence, prob) {
    seq_chars <- strsplit(sequence, "")[[1]]
    lowercase_positions <- sample(1:length(seq_chars),
                                  size = round(length(seq_chars) * prob))
    seq_chars[lowercase_positions] <- tolower(seq_chars[lowercase_positions])
    return(paste(seq_chars, collapse = ""))
  }

  # Apply to all sequences
  data$AlleleSequence <- sapply(data$AlleleSequence,
                                function(x) convert_sequence(x, lowercase_prob))
  return(data)
}

# Simple function to add random IUPAC codes to sequences
add_random_iupac <- function(data, iupac_prob = 0.05) {
  # Function to convert random bases to IUPAC codes in a single sequence
  convert_sequence <- function(sequence, prob) {
    seq_chars <- strsplit(sequence, "")[[1]]

    # Select random positions to convert to IUPAC
    n_positions <- round(length(seq_chars) * prob)
    if (n_positions > 0) {
      iupac_positions <- sample(1:length(seq_chars), size = n_positions)

      # Available IUPAC codes (Y, R, and assuming you meant N instead of I)
      iupac_codes <- c("Y", "R", "N")  # Y=C/T, R=A/G, N=any base

      # Replace selected positions with random IUPAC codes
      seq_chars[iupac_positions] <- sample(iupac_codes, size = n_positions, replace = TRUE)
    }

    return(paste(seq_chars, collapse = ""))
  }

  # Apply to all sequences
  data$AlleleSequence <- sapply(data$AlleleSequence,
                                function(x) convert_sequence(x, iupac_prob))
  return(data)
}


# Main execution
# Update these paths to your actual files
fasta_file_path <- "~/Documents/github/BIGapp-PanelHub/alfalfa/alfalfa_allele_db_v001.fa"
snp_file_path <- "~/Documents/github/BIGapp-PanelHub/alfalfa/20201030-BI-Alfalfa_SNPs_DArTag-probe-design_snpID_lut.csv"  # Your SNP information file

# Check if files exist
if (!file.exists(fasta_file_path)) {
  stop("FASTA file not found. Please check the file path: ", fasta_file_path)
}
if (!file.exists(snp_file_path)) {
  stop("SNP file not found. Please check the file path: ", snp_file_path)
}

# Generate the data
cat("Starting simulation with real SNP data and realistic count distributions...\n")
set.seed(123)  # For reproducibility

# Specify how many CloneIDs you want (change this number as needed)
n_clones_wanted <- 300  # Change this to any number you want
n_samples_wanted <- 30  # Change this if you want different number of samples

simulated_data <- generate_allele_data(
  fasta_file = fasta_file_path,
  snp_file = snp_file_path,
  n_clones = n_clones_wanted,
  n_samples = n_samples_wanted
)

# Write to file
output_file <- "alfalfa_madc.csv"
write.csv(simulated_data, output_file, row.names = FALSE)

# WRong ID

simulated_data1 <- simulated_data
simulated_data1$CloneID <- gsub("[.]1", "",simulated_data$CloneID)

output_file <- "alfalfa_madc_wrongID.csv"
write.csv(simulated_data1, output_file, row.names = FALSE)

cat("Simulation complete!\n")
cat("Generated file:", output_file, "\n")
cat("Total rows:", nrow(simulated_data), "\n")
cat("Total columns:", ncol(simulated_data), "\n")
cat("Unique CloneIDs:", length(unique(simulated_data$CloneID)), "\n")

# Display first few rows
cat("\nFirst few rows preview:\n")
print(head(simulated_data[, 1:min(8, ncol(simulated_data))], 10))

# Summary statistics
cat("\nAlleleID types summary:\n")
allele_types <- gsub(".*\\|(.+?)(_.*)?$", "\\1", simulated_data$AlleleID)
print(table(allele_types))

# Check zero distribution
sample_cols <- grep("Sample_", names(simulated_data), value = TRUE)
zero_counts <- sapply(simulated_data[sample_cols], function(x) sum(x == 0))
cat("\nZero count distribution across samples:\n")
cat("Mean zeros per sample:", round(mean(zero_counts), 1), "\n")
cat("Range of zeros per sample:", min(zero_counts), "to", max(zero_counts), "\n")
cat("Total zeros in dataset:", sum(zero_counts), "out of", nrow(simulated_data) * length(sample_cols), "observations\n")
cat("Percentage of zeros:", round(sum(zero_counts) / (nrow(simulated_data) * length(sample_cols)) * 100, 1), "%\n")


# Creating examples for the other scenarios
## Lower case

# Add random lowercase
simulated_datal <- add_random_lowercase(simulated_data, lowercase_prob = 0.15)
write.csv(simulated_datal, "alfalfa_lowercase.csv", row.names = FALSE)

# Check the result
head(simulated_datal$AlleleSequence)

## IUPAC codes

# Add random IUPAC codes
simulated_data <- add_random_iupac(simulated_data, iupac_prob = 0.05)

# Check the result
head(simulated_data$AlleleSequence)
write.csv(simulated_data, "alfalfa_IUPAC.csv", row.names = FALSE)

## With indel and bad ChromPOS

fasta_file_path = "~/Documents/github/BIGapp-PanelHub/potato/potato_allele_db_v001.fa"
snp_file_path = "~/Documents/github/BIGapp-PanelHub/potato/potato_dartag_v2_3915markers_rm7dupTags_6traitMarkers_rm1dup_snpID_lut.csv"

specific_clones <- c("chr05_004488021", "chr05_004488015")  # Your specific CloneIDs

set.seed(1213)
simulated_data <- generate_allele_data(
  fasta_file = fasta_file_path,
  snp_file = snp_file_path,
  n_clones = 300,
  n_samples = 30,
  force_include = specific_clones  # These will always be included
)

# Write to file
output_file <- "potato_indel_madc.csv"
write.csv(simulated_data, output_file, row.names = FALSE)

## With non-standard CloneID
mk_info <- read.csv("~/Documents/github/BIGapp-PanelHub/potato/potato_dartag_v2_3915markers_rm7dupTags_6traitMarkers_rm1dup_snpID_lut.csv")
simulated_data$CloneID <- mk_info$Panel_markerID[match(simulated_data$CloneID, mk_info$BI_markerID)]
write.csv(simulated_data, "potato_indel_ChromPosFALSE.csv", row.names = FALSE)

## with just some non-standard CloneID and more indels
specific_clones <- c("chr05_004488021", "chr05_004488015",
                     "M6_chr10_48867893_000000440",
                     "M6_chr10_48867893_000000225",
                     "C88_C10H2G055580_399_000000151")  # Your specific CloneIDs

set.seed(1231)
simulated_data <- generate_allele_data(
  fasta_file = fasta_file_path,
  snp_file = snp_file_path,
  n_clones = 300,
  n_samples = 30,
  force_include = specific_clones  # These will always be included
)

# Write to file
output_file <- "potato_more_indels_madc_ChromPosFALSE.csv"
write.csv(simulated_data, output_file, row.names = FALSE)

## Lower case

# Add random lowercase
simulated_data1 <- add_random_lowercase(simulated_data, lowercase_prob = 0.15)
write.csv(simulated_data1, "potato_indel_lowercase.csv", row.names = FALSE)

# Check the result
head(simulated_data$AlleleSequence)

## IUPAC codes

# Add random IUPAC codes
simulated_data <- add_random_iupac(simulated_data, iupac_prob = 0.05)

# Check the result
head(simulated_data$AlleleSequence)
write.csv(simulated_data, "potato_indel_IUPAC.csv", row.names = FALSE)

