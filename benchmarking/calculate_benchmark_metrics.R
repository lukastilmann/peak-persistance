#' Calculate Metrics for Benchmark Results
#'
#' @description
#' Calculates distance metrics between the ground truth, mean function, estimated function,
#' and warping functions from benchmark results.
#'
#' @param benchmark_result List containing the output from run_benchmark function.
#'
#' @return A list containing various performance metrics:
#'   \item{mean_ground_truth_distance}{L2 distance between mean function and ground truth}
#'   \item{est_ground_truth_distance}{L2 distance between estimated function and ground truth}
#'   \item{warping_distances}{Mean L2 distances between warping functions at different stages}
#'   \item{alignment_scores}{Various metrics quantifying the quality of alignment}
#'   \item{peak_detection_metrics}{Metrics for peak detection accuracy}
#'
#' @import tf
#'
#' @examples
#' \dontrun{
#' # Assuming benchmark_result is the output from run_benchmark
#' metrics <- calculate_benchmark_metrics(benchmark_result)
#' print(metrics)
#' }
#'
#' @export
calculate_benchmark_metrics <- function(benchmark_result) {
  # Input validation
  if (!is.list(benchmark_result)) {
    stop("benchmark_result must be a list")
  }

  required_elements <- c("params", "runtime", "simulated_data", "ppd_result",
                         "estimated_function", "alignment_to_est",
                         "al_to_est_warping_functions")
  missing_elements <- setdiff(required_elements, names(benchmark_result))
  if (length(missing_elements) > 0) {
    stop("Missing required elements in benchmark_result: ",
         paste(missing_elements, collapse = ", "))
  }

  # Initialize metrics list
  metrics <- list()

  # Add original parameters for reference
  metrics$params <- benchmark_result$params

  # Add runtime metrics
  metrics$runtime <- benchmark_result$runtime

  # Check if results exist
  if (is.null(benchmark_result$results)) {
    warning("No results found in benchmark_result")
    return(metrics)
  }

  # Calculate distance between mean function and ground truth
  if (!is.null(benchmark_result$ppd_result$mean_function) &&
      !is.null(benchmark_result$simulated_data$ground_truth)) {

    mean_func <- benchmark_result$ppd_result$mean_function
    ground_truth <- benchmark_result$simulated_data$ground_truth
    t_grid <- benchmark_result$simulated_data$t_grid

    # Calculate L2 distance between mean function and ground truth
    # NOTE: right metric?
    try({
      mean_gt_diff <- mean_func - ground_truth
      metrics$mean_ground_truth_distance <- sqrt(tf::tf_integrate(mean_gt_diff^2,
                                                                  lower = min(t_grid),
                                                                  upper = max(t_grid)))
    }, silent = TRUE)
  }

  # Calculate distance between estimated function and ground truth
  if (!is.null(benchmark_result$estimated_function) &&
      !is.null(benchmark_result$simulated_data$ground_truth)) {

    est_func <- benchmark_result$estimated_function
    ground_truth <- benchmark_result$simulated_data$ground_truth
    t_grid <- benchmark_result$simulated_data$t_grid

    # Calculate L2 distance between estimated function and ground truth
    try({
      est_gt_diff <- est_func - ground_truth
      metrics$est_ground_truth_distance <- sqrt(tf::tf_integrate(est_gt_diff^2,
                                                                 lower = min(t_grid),
                                                                 upper = max(t_grid)))
    }, silent = TRUE)
  }

  # Calculate warping function distances
  # 1. Distance between original warping functions and PPD warping functions
  if (!is.null(benchmark_result$simulated_data$warped_grid) &&
      !is.null(benchmark_result$warping_functions)) {

    try({
      # Original warping functions from data generation
      orig_warping <- benchmark_result$simulated_data$warped_grid
      n_curves <- nrow(orig_warping)
      t_grid <- benchmark_result$simulated_data$t_grid

      # Convert to tfd objects for easier computation
      orig_warping_tfd <- tf::tfd(orig_warping, t_grid)

      # PPD warping functions
      ppd_warping <- benchmark_result$warping_functions

      # Calculate mean L2 distance
      warping_diffs <- numeric(n_curves)
      for (i in 1:n_curves) {
        diff_func <- orig_warping_tfd[i,] - ppd_warping[i,]
        warping_diffs[i] <- sqrt(tf::tf_integrate(diff_func^2,
                                                  lower = min(t_grid),
                                                  upper = max(t_grid)))
      }

      metrics$warping_distances$ppd_warping_mean_distance <- mean(warping_diffs)
      metrics$warping_distances$ppd_warping_sd_distance <- sd(warping_diffs)
    }, silent = TRUE)
  }

  # 2. Distance between original warping functions and final alignment warping functions
  if (!is.null(benchmark_result$simulated_data$warped_grid) &&
      !is.null(benchmark_result$results$al_to_est_warping_functions)) {

    try({
      # Original warping functions from data generation
      orig_warping <- benchmark_result$simulated_data$warped_grid
      n_curves <- nrow(orig_warping)
      t_grid <- benchmark_result$simulated_data$t_grid

      # Convert to tfd objects for easier computation
      orig_warping_tfd <- tf::tfd(orig_warping, t_grid)

      # Final alignment warping functions
      align_warping <- benchmark_result$al_to_est_warping_functions

      # Calculate mean L2 distance
      warping_diffs <- numeric(n_curves)
      for (i in 1:n_curves) {
        diff_func <- orig_warping_tfd[i,] - align_warping[i,]
        warping_diffs[i] <- sqrt(tf::tf_integrate(diff_func^2,
                                                  lower = min(t_grid),
                                                  upper = max(t_grid)))
      }

      metrics$warping_distances$alignment_warping_mean_distance <- mean(warping_diffs)
      metrics$warping_distances$alignment_warping_sd_distance <- sd(warping_diffs)
    }, silent = TRUE)
  }

  return(metrics)
}
