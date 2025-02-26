#' Track Peak Succession Between Curves
#'
#' @description
#' Identifies which peaks in a second curve correspond to peaks in a first curve,
#' maintaining consistent labeling across different lambda values.
#'
#' @param peaks_1 Numeric vector of peak locations in the first curve
#' @param valleys_1 Numeric vector of valley locations in the first curve
#' @param peaks_2 Numeric vector of peak locations in the second curve
#' @param labels_1 Numeric vector of labels for peaks in the first curve
#' @param label_max Numeric value indicating the current maximum label used
#' @param t_grid Numeric vector representing the time grid for the functional data
#'
#' @return A list containing:
#' \itemize{
#'   \item labels: Numeric vector of labels for peaks in the second curve
#'   \item label_max: Updated maximum label value
#' }
#'
#' @details
#' The function tracks peaks across different curves by determining if peaks in the
#' second curve fall within the range defined by valleys surrounding peaks in the first curve.
peak_successor <- function(peaks_1, valleys_1, peaks_2, labels_1,
                           label_max, t_grid) {
  # Input validation
  if (missing(peaks_1)) {
    stop("'peaks_1' is required")
  }
  if (missing(valleys_1)) {
    stop("'valleys_1' is required")
  }
  if (missing(peaks_2)) {
    stop("'peaks_2' is required")
  }
  if (missing(labels_1)) {
    stop("'labels_1' is required")
  }
  if (missing(label_max)) {
    stop("'label_max' is required")
  }
  if (missing(t_grid)) {
    stop("'t_grid' is required")
  }

  # Type checking
  if (!is.numeric(peaks_1)) {
    stop("'peaks_1' must be a numeric vector")
  }
  if (!is.numeric(valleys_1)) {
    stop("'valleys_1' must be a numeric vector")
  }
  if (!is.numeric(peaks_2)) {
    stop("'peaks_2' must be a numeric vector")
  }
  if (!is.numeric(labels_1)) {
    stop("'labels_1' must be a numeric vector")
  }
  if (!is.numeric(label_max) || length(label_max) != 1) {
    stop("'label_max' must be a single numeric value")
  }
  if (!is.numeric(t_grid)) {
    stop("'t_grid' must be a numeric vector")
  }

  # Consistency checks
  if (length(peaks_1) != length(labels_1)) {
    stop("'peaks_1' and 'labels_1' must have the same length")
  }
  if (length(t_grid) < 2) {
    stop("'t_grid' must have at least two points")
  }
  if (!all(is.finite(t_grid))) {
    stop("'t_grid' must contain only finite values")
  }

  # Helper function to compute ranges for peaks
  compute_peak_ranges <- function(peaks, valleys) {
    # Ensure valleys include endpoints
    all_valleys <- sort(unique(c(t_grid[1], valleys, t_grid[length(t_grid)])))
    # For each peak, find its containing range
    ranges <- lapply(peaks, function(peak) {
      left_valley <- max(all_valleys[all_valleys < peak])
      right_valley <- min(all_valleys[all_valleys > peak])
      c(left_valley, right_valley)
    })
    # Convert list to matrix
    do.call(rbind, ranges)
  }

  # Helper function to assign labels to peaks based on ranges
  assign_labels_to_peaks <- function(peaks_2, ranges, peaks_1, labels_1) {
    labels <- numeric(length(peaks_2))
    for(i in seq_along(peaks_2)) {
      # Find which ranges contain this peak
      in_range <- peaks_2[i] >= ranges[,1] & peaks_2[i] <= ranges[,2]
      matching_ranges <- which(in_range)
      if(length(matching_ranges) > 0) {
        if(length(matching_ranges) > 1) {
          # If multiple ranges match, choose the closest peak
          distances <- abs(peaks_1[matching_ranges] - peaks_2[i])
          closest_idx <- matching_ranges[which.min(distances)]
          labels[i] <- labels_1[closest_idx]
        } else {
          labels[i] <- labels_1[matching_ranges]
        }
      }
    }
    labels
  }

  # Helper function to resolve overlapping labels
  resolve_overlapping_labels <- function(labels, peaks_1, peaks_2, labels_1) {
    unique_labels <- unique(labels[labels > 0])
    for(label in unique_labels) {
      duplicates <- which(labels == label)
      if(length(duplicates) > 1) {
        # Keep the closest peak and reset others
        orig_peak <- peaks_1[which(labels_1 == label)]
        distances <- abs(orig_peak - peaks_2[duplicates])
        keep_idx <- duplicates[which.min(distances)]
        remove_idx <- duplicates[duplicates != keep_idx]
        labels[remove_idx] <- 0
      }
    }
    labels
  }

  # Main function logic
  if(length(peaks_1) == 0) {
    # If no peaks in first function, assign new labels to all peaks in second
    labels_2 <- label_max + seq_along(peaks_2)
    return(list(
      labels = labels_2,
      label_max = label_max + length(peaks_2)
    ))
  }

  # Compute ranges for peaks in first function
  ranges <- compute_peak_ranges(peaks_1, valleys_1)

  # Assign initial labels to peaks in second function
  labels_2 <- assign_labels_to_peaks(peaks_2, ranges, peaks_1, labels_1)

  # Resolve any overlapping labels
  labels_2 <- resolve_overlapping_labels(labels_2, peaks_1, peaks_2, labels_1)

  # Assign new labels to unmatched peaks
  unmatched <- which(labels_2 == 0)
  if(length(unmatched) > 0) {
    labels_2[unmatched] <- label_max + seq_along(unmatched)
    label_max <- label_max + length(unmatched)
  }

  # Return both new labels and updated label_max
  return(list(
    labels = labels_2,
    label_max = label_max
  ))
}
