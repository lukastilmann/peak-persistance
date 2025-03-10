#' Extract Lambda Value from Column Name
#'
#' @param col_name String column name containing a lambda value
#' @return Numeric lambda value
get_lambda <- function(col_name) {
  # Check inputs
  if (missing(col_name)) {
    stop("'col_name' is required")
  }
  if (!(is.character(col_name) || is.factor(col_name))) {
    stop("'col_name' must be a character string")
  }

  # Extract lambda value
  as.numeric(sub("aligned_lambda_", "", col_name))
}


#' Get Significant Peaks
#'
#' @param curvatures Numeric vector of peak curvatures
#' @param labels Numeric vector of peak labels
#' @param tau Significance threshold for curvatures
#'
#' @return Numeric vector of significant peak labels
get_significant_peaks <- function(curvatures, labels, tau) {
  # Input validation
  if (missing(curvatures)) {
    stop("'curvatures' is required")
  }
  if (missing(labels)) {
    stop("'labels' is required")
  }
  if (missing(tau)) {
    stop("'tau' is required")
  }

  if (!is.numeric(curvatures)) {
    stop("'curvatures' must be numeric")
  }
  if (!is.numeric(labels)) {
    stop("'labels' is required")
  }
  if (!is.numeric(tau) || length(tau) != 1) {
    stop("'tau' must be a single numeric value")
  }

  if (length(curvatures) != length(labels)) {
    stop("'curvatures' and 'labels' must have the same length")
  }

  sig_curvatures <- unlist(curvatures > tau)
  sig_labels <- labels[sig_curvatures]
  sig_labels
}

#' Find Optimal Lambda Index
#'
#' @param significant_peaks List of significant peak labels for each lambda
#' @param persistent_labels Numeric vector of persistent peak labels
#'
#' @return Index of optimal lambda value
find_optimal_lambda <- function(significant_peaks, persistent_labels) {
  # Input validation
  if (missing(significant_peaks)) {
    stop("'significant_peaks' is required")
  }
  if (missing(persistent_labels)) {
    stop("'persistent_labels' is required")
  }

  if (!is.list(significant_peaks)) {
    stop("'significant_peaks' must be a list")
  }
  if (!is.numeric(persistent_labels)) {
    stop("'persistent_labels' must be numeric")
  }

  if (length(significant_peaks) == 0) {
    stop("'significant_peaks' cannot be empty")
  }

  # Flatten all peaks
  all_peaks <- unlist(significant_peaks)
  if (length(all_peaks) == 0) {
    return(1)  # Default to first lambda if no significant peaks
  }

  # Check if any persistent labels are outside the range
  if (length(persistent_labels) > 0 && max(persistent_labels) > max(all_peaks)) {
    stop("'persistent_labels' contains values outside the range of peak labels")
  }

  # Convert significant peaks to matrix form
  n_lams <- length(significant_peaks)
  ref_row <- numeric(max(all_peaks))
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

#' Get Peak Intervals
#'
#' @param significant_peaks List of significant peak labels for each lambda
#'
#' @return List where names are peak labels and values are lambda values where the peak is significant
get_peak_intervals <- function(significant_peaks) {
  # Input validation
  if (missing(significant_peaks)) {
    stop("'significant_peaks' is required")
  }

  if (!is.list(significant_peaks)) {
    stop("'significant_peaks' must be a list")
  }

  if (length(significant_peaks) == 0) {
    stop("'significant_peaks' cannot be empty")
  }

  if (is.null(names(significant_peaks))) {
    stop("'significant_peaks' must be a named list")
  }

  # Check if get_lambda function exists
  if (!exists("get_lambda")) {
    stop("Function 'get_lambda' is required but not available")
  }

  # Create a map of all unique labels and their appearances
  all_labels <- unique(unlist(significant_peaks))

  # Handle empty result
  if (length(all_labels) == 0) {
    return(list())
  }

  # Get lambda values from names
  lambda_values <- tryCatch({
    unlist(lapply(names(significant_peaks), function(name) {
      get_lambda(name)
    }))
  }, error = function(e) {
    stop("Error extracting lambda values: ", e$message)
  })

  # Initialize list to store intervals for each label
  label_intervals <- list()

  # For each label, find all lambda values where it appears
  for (label in all_labels) {
    # Find in which lambda indices this label appears
    appearances <- which(sapply(significant_peaks, function(peaks) label %in% peaks))
    if (length(appearances) > 0) {
      label_intervals[[as.character(label)]] <- lambda_values[appearances]
    }
  }

  return(label_intervals)
}


#' Format Lambda Values for Display
#'
#' @description Formats a vector of lambda values for better readability,
#'   using scientific notation for very small values and fixed notation for others.
#'
#' @param lambda_values Numeric vector of lambda values to format
#'
#' @return Character vector of formatted lambda values
#'
#' @details Values less than 0.0001 are formatted using scientific notation with 4 digits,
#'   while larger values use fixed-point notation with 4 decimal places.
#'
format_lambda <- function(lambda_values) {
  result <- character(length(lambda_values))
  for (i in seq_along(lambda_values)) {
    # Check if the value is very small (< 0.0001)
    if (abs(lambda_values[i]) < 0.0001) {
      # Use scientific notation for very small values
      result[i] <- formatC(lambda_values[i], format = "e", digits = 4)
    } else {
      # Use fixed notation for other values
      result[i] <- formatC(lambda_values[i], format = "f", digits = 4)
    }
  }
  return(result)
}
