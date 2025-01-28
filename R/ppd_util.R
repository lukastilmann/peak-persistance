
# Finds which peaks are significant (curvature above threshold)
get_significant_peaks <- function(curvatures, labels, tau){
  sig_curvatures <- unlist(curvatures > tau)
  sig_labels <- labels[sig_curvatures]
  sig_labels
}


# Finds which peaks are significant for a large enough interval
find_persistent_peaks <- function(significant_peaks, ratio_threshold) {
  # Create a map of all unique labels and their appearances
  all_labels <- unique(unlist(significant_peaks))

  # Get lambda values from names
  lambda_values <- unlist(lapply(names(significant_peaks), function(name){
    get_lambda(name)
  }))

  # Initialize list to store intervals for each label
  label_intervals <- list()

  # For each label, find all lambda values where it appears
  for(label in all_labels) {
    # Find in which lambda indices this label appears
    appearances <- which(sapply(significant_peaks, function(peaks) label %in% peaks))

    if(length(appearances) > 0) {
      label_intervals[[as.character(label)]] <- lambda_values[appearances]
    }
  }

  # Find the label with the most appearances
  max_appearances <- max(sapply(label_intervals, length))

  # Find persistent peaks (those appearing in >= ratio_threshold * max_appearances lambdas)
  persistent_labels <- names(which(sapply(label_intervals, length) >= ratio_threshold * max_appearances))
  persistent_labels <- as.numeric(persistent_labels)

  return(list(
    persistent_labels = persistent_labels,
    intervals = label_intervals
  ))
}


# Function to find optimal lambda index
find_optimal_lambda <- function(significant_peaks, persistent_labels) {
  # Convert significant peaks to matrix form
  n_lams <- length(significant_peaks)
  ref_row <- numeric(max(unlist(significant_peaks)))
  ref_row[persistent_labels] <- 1

  # Calculate Hamming distances
  hamming_distances <- sapply(significant_peaks, function(peaks) {
    comp <- numeric(length(ref_row))
    comp[peaks] <- 1
    sum(comp != ref_row)
  })

  # Find exact matches
  exact_matches <- which(hamming_distances == 0)

  if (length(exact_matches) > 0) {
    return(min(exact_matches))
  } else {
    min_distance <- min(hamming_distances)
    return(min(which(hamming_distances == min_distance)))
  }
}
