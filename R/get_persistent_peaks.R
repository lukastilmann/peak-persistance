#' Find Persistent Peaks Using a Threshold Method
#'
#' @param label_intervals List of peak label intervals
#' @param ratio_threshold Numeric threshold between 0 and 1 for determining persistence
#'
#' @return Numeric vector of persistent peak labels
get_persistent_peaks_threshold <- function(label_intervals, ratio_threshold) {
  # Input validation
  if (missing(label_intervals)) {
    stop("'label_intervals' is required")
  }
  if (!is.list(label_intervals)) {
    stop("'label_intervals' must be a list")
  }
  if (length(label_intervals) == 0) {
    return(numeric(0))
  }
  if (missing(ratio_threshold)) {
    stop("'ratio_threshold' is required")
  }
  if (!is.numeric(ratio_threshold) || length(ratio_threshold) != 1) {
    stop("'ratio_threshold' must be a single numeric value")
  }
  if (ratio_threshold < 0 || ratio_threshold > 1) {
    stop("'ratio_threshold' must be between 0 and 1")
  }

  # Find the label with the most appearances
  max_appearances <- max(sapply(label_intervals, length))

  # Find persistent peaks (those appearing in >= ratio_threshold * max_appearances lambdas)
  persistent_labels <- names(which(sapply(label_intervals, length) >= ratio_threshold * max_appearances))
  persistent_labels <- as.numeric(persistent_labels)

  return(persistent_labels)
}


#' Find Persistent Peaks Using Clustering Method
#'
#' @param peak_lists List of significant peaks for each lambda value
#'
#' @return Numeric vector of persistent peak labels
get_persistent_peaks_clustering <- function(peak_lists) {
  # Input validation
  if (missing(peak_lists)) {
    stop("'peak_lists' is required")
  }
  if (!is.list(peak_lists)) {
    stop("'peak_lists' must be a list")
  }
  if (length(peak_lists) == 0) {
    return(numeric(0))
  }

  # Check if each element of the list is a numeric vector
  if (!all(sapply(peak_lists, is.numeric))) {
    stop("All elements in 'peak_lists' must be numeric vectors")
  }

  # Check if peak_lists is empty or contains only empty vectors
  if (all(sapply(peak_lists, length) == 0)) {
    return(numeric(0))
  }

  # Create indicator matrix
  all_peaks <- unlist(peak_lists)
  n_lambdas <- length(peak_lists)
  max_peak <- max(all_peaks)
  row_idx <- rep(seq_len(n_lambdas), lengths(peak_lists))
  col_idx <- all_peaks
  indices <- matrix(c(row_idx, col_idx), ncol = 2)
  indicator_matrix <- matrix(0, nrow = n_lambdas, ncol = max_peak)
  indicator_matrix[indices] <- 1

  # Count the number of ones (occurrences) for each peak (ignore NAs)
  occurence_counts <- colSums(indicator_matrix, na.rm = TRUE)

  # Add dummy peak with no occurences
  occ_data <- c(occurence_counts, 0)
  if (all(occ_data == 0) || length(unique(occurence_counts)) == 1) {
    return(numeric(0))
  }

  # Compute pairwise distances and perform hierarchical clustering
  distances <- stats::dist(occ_data)
  hierarchical_clustering <- stats::hclust(distances, method = "ward.D2")

  # Cluster the occ_data into 2 clusters
  cluster_assignments <- stats::cutree(hierarchical_clustering, k = 2)

  # Get reference cluster from the last element (representing peak with no occurences)
  reference_cluster <- cluster_assignments[length(cluster_assignments)]
  cluster_assignments <- cluster_assignments[-length(cluster_assignments)]

  # Identify indices where cluster assignments differ from reference
  persistent_peaks <- which(cluster_assignments != reference_cluster)

  return(persistent_peaks)
}
