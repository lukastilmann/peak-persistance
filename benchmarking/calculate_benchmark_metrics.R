#' Calculate Metrics for Benchmark Results
#'
#' @description
#' Calculates comprehensive distance metrics between the ground truth, mean function, estimated function,
#' and warping functions from benchmark results. These metrics quantify alignment quality and functional
#' approximation accuracy across different stages of the analysis.
#'
#' @param benchmark_result List containing the output from run_benchmark function. Must include components
#'   such as params, runtime, simulated_data, ppd_result, estimated_function, alignment_to_est, and
#'   al_to_est_warping_functions.
#'
#' @return A list containing various performance metrics:
#'   \item{params}{Original parameters from the benchmark}
#'   \item{runtime}{Runtime measurements from the benchmark process}
#'   \item{mean_func_ground_truth_dist}{L2 distance between mean function and ground truth}
#'   \item{mean_func_amplitude_dist}{Amplitude distance between aligned mean function and ground truth}
#'   \item{mean_func_phase_dist}{Phase distance between mean function and ground truth}
#'   \item{mean_func_gt_warping_dist}{L2 distance between mean function warping and identity warping}
#'   \item{sce_func_ground_truth_dist}{L2 distance between estimated function and ground truth}
#'   \item{sce_func_amplitude_dist}{Amplitude distance between aligned estimated function and ground truth}
#'   \item{sce_func_phase_dist}{Phase distance between estimated function and ground truth}
#'   \item{sce_func_gt_warping_dist}{L2 distance between estimated function warping and identity warping}
#'   \item{warping_distances}{Mean and standard deviation of L2 distances between warping functions at different stages}
#'     \item{ppd_warping_mean_distance}{Mean L2 distance between original and PPD warping functions}
#'     \item{ppd_warping_sd_distance}{Standard deviation of L2 distances between original and PPD warping functions}
#'     \item{alignment_warping_mean_distance}{Mean L2 distance between original and final alignment warping functions}
#'     \item{alignment_warping_sd_distance}{Standard deviation of L2 distances between original and final alignment warping functions}}
#'   \item{num_peaks}{Number of peaks detected in the estimated function}
#'   \item{num_peaks_error}{Difference between detected number of peaks and true number of peaks}
#'
#' @import tf
#'
#' @examples
#' \dontrun{
#' # Assuming benchmark_result is the output from run_benchmark
#' metrics <- calculate_benchmark_metrics(benchmark_result)
#' print(metrics)
#'
#' # View specific metrics
#' print(metrics$mean_ground_truth_distance)
#' print(metrics$warping_distances)
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

  # Accessing results needed for multiple computations
  if (!is.null(benchmark_result$simulated_data$ground_truth) &&
      !is.null(benchmark_result$simulated_data$t_grid)){

    ground_truth <- benchmark_result$simulated_data$ground_truth
    t_grid <- benchmark_result$simulated_data$t_grid
    try({
      identity_warp <- tf::tfd(t_grid, t_grid)
    })
  }

  # Calculate distance between mean function and ground truth
  if (!is.null(benchmark_result$ppd_result$mean_function)) {

    mean_func <- benchmark_result$ppd_result$mean_function

    # Calculate L2 distance between mean function and ground truth
    try({
      mean_gt_diff <- mean_func - ground_truth
      metrics$mean_func_ground_truth_dist <- sqrt(tf::tf_integrate(mean_gt_diff^2,
                                                                  lower = min(t_grid),
                                                                  upper = max(t_grid)))
    }, silent = TRUE)

    # Align mean function to ground truth function
    try({
      res_align <- align_pairwise(ground_truth, mean_func, t_grid)
      q_1 <- res_align$q_1
      q_2_aligned <- res_align$q_2_aligned
      warping_function <- res_align$warping_function

      amplitude_dist_mean <- sqrt(tf::tf_integrate((q_1 - q_2_aligned)^2,
                           lower = tf::tf_domain(q_1)[1], upper = tf_domain(q_1)[2]))

      # Convert to basis representation to compute derivative
      warping_function_b <- tf::tfb(warping_function, arg = t_grid, basis = "spline")
      # Calculating phase distance metric
      warping_function_deriv <- tf::tfd(tf::tf_derive(warping_function_b), t_grid)
      warping_function_deriv_sqrt <- sqrt(warping_function_deriv)
      warp_deriv_integral <- tf::tf_integrate(warping_function_deriv_sqrt,
                                              lower = tf::tf_domain(warping_function_deriv_sqrt)[1],
                                              upper = tf::tf_domain(warping_function_deriv_sqrt)[2])
      phase_dist_mean <- acos(min(warp_deriv_integral, 1))


      # Distance of warping function from identity warping function
      mfn_gt_diff <- warping_function - identity_warp
      mean_func_gt_warping_distance <- sqrt(tf::tf_integrate(mfn_gt_diff^2,
                                                             lower = tf::tf_domain(mfn_gt_diff)[1],
                                                             upper = tf::tf_domain(mfn_gt_diff)[2]))

      # Save metrics
      metrics$mean_func_amplitude_dist <- amplitude_dist_mean
      metrics$mean_func_phase_dist <- phase_dist_mean
      metrics$mean_func_gt_warping_dist <- mean_func_gt_warping_distance

    }, silent = TRUE)
  }

  # Calculate distance between estimated function and ground truth
  if (!is.null(benchmark_result$estimated_function)) {

    est_func <- benchmark_result$estimated_function

    # Calculate L2 distance between estimated function and ground truth
    try({
      est_gt_diff <- est_func - ground_truth
      metrics$sce_func_ground_truth_dist <- sqrt(tf::tf_integrate(est_gt_diff^2,
                                                                 lower = min(t_grid),
                                                                 upper = max(t_grid)))

      # TODO: alignment and phase amplitude distance
    }, silent = TRUE)

    # Align shape constrained estimate function to ground truth function
    try({
      res_align_sce <- align_pairwise(ground_truth, est_func, t_grid)
      q_1 <- res_align_sce$q_1
      q_2_aligned <- res_align_sce$q_2_aligned
      warping_function <- res_align_sce$warping_function

      amplitude_dist_sce <- sqrt(tf::tf_integrate((q_1 - q_2_aligned)^2,
                                                   lower = tf::tf_domain(q_1)[1],
                                                   upper = tf::tf_domain(q_1)[2]))

      # Convert to basis representation to compute derivative
      warping_function_b <- tf::tfb(warping_function, arg = t_grid, basis = "spline")
      # Calculating phase distance metric
      warping_function_deriv <- tf::tfd(tf::tf_derive(warping_function_b), t_grid)
      warping_function_deriv_sqrt <- sqrt(warping_function_deriv)
      warp_deriv_integral <- tf::tf_integrate(warping_function_deriv_sqrt,
                                              lower = tf::tf_domain(warping_function_deriv_sqrt)[1],
                                              upper = tf::tf_domain(warping_function_deriv_sqrt)[2])
      phase_dist_sce <- acos(min(warp_deriv_integral, 1))

      # Distance of warping function from identity warping function
      sce_gt_diff <- warping_function - identity_warp
      sce_func_gt_warping_distance <- sqrt(tf::tf_integrate(sce_gt_diff^2,
                                                             lower = tf::tf_domain(sce_gt_diff)[1],
                                                             upper = tf::tf_domain(sce_gt_diff)[2]))

      metrics$sce_func_amplitude_dist <- amplitude_dist_sce
      metrics$sce_func_phase_dist <- phase_dist_sce
      metrics$sce_func_gt_warping_dist <- sce_func_gt_warping_distance

    }, silent = TRUE)
  }

  # Calculate warping function distances
  # 1. Distance between original warping functions and PPD warping functions
  if (!is.null(benchmark_result$simulated_data$warped_grid) &&
      !is.null(benchmark_result$ppd_result$warping_functions)) {

    try({
      # Original warping functions from data generation
      orig_warping <- benchmark_result$simulated_data$warped_grid
      n_curves <- nrow(orig_warping)
      t_grid <- benchmark_result$simulated_data$t_grid

      # Convert to tfd objects for easier computation
      orig_warping_tfd <- tf::tfd(orig_warping, t_grid)

      # PPD warping functions
      ppd_warping <- benchmark_result$ppd_result$warping_functions

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
      !is.null(benchmark_result$al_to_est_warping_functions)) {

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

  # Error in number of peaks
  if (!is.null(benchmark_result$ppd_result$num_peaks) &&
      !is.null(benchmark_result$true_num_peaks)) {

    try({
      metrics$num_peaks <- benchmark_result$ppd_result$num_peaks
      metrics$num_peaks_error <- benchmark_result$ppd_result$num_peaks - benchmark_result$true_num_peaks
    }, silent = TRUE)
  }

  return(metrics)
}
