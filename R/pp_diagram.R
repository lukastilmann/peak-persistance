#' Create Peak Persistence Diagram for Aligned Functions
#'
#' @description
#' This function creates barchart and surface diagrams tracking peaks of aligned functions
#' for different lambda values. It identifies significant and persistent peaks across
#' multiple levels of alignment regularization.
#'
#' @param curves An object containing functional data to be analyzed, compatible with the tfb function.
#' @param t_grid Numeric vector representing the time grid for the functional data.
#' @param parallel Logical indicating whether to use parallel processing.
#' @param lambda_search_start Numeric, starting value for searching maximum lambda (default: 2).
#' @param lambda_search_min_bound Numeric, minimum bound for lambda search (default: 0.01).
#' @param lambda_search_threshold Numeric, convergence threshold for lambda search (default: 1e-3).
#' @param max_lambda_search_steps Integer, maximum number of steps for lambda search (default: 20).
#' @param max_iter Integer, maximum number of iterations for align_functions (default: 20).
#' @param penalty Character, penalty used in function alignment.
#' @param n_lambda Integer, number of lambda values to evaluate (default: 10).
#' @param sig_threshold Numeric, threshold for peak significance (default: 0.03).
#' @param pers_threshold Numeric, threshold for peak persistence when using 'threshold' method (default: 0.28).
#' @param pers_method Character, method for determining persistent peaks: "clustering" or "threshold" (default: "clustering").
#'
#' @return A list containing:
#' \itemize{
#'   \item bc: Persistence diagram in barchart form
#'   \item surface: Peak persistence surface plot
#'   \item significant_peaks: Vector of significant peaks at optimal lambda
#'   \item persistent_peaks: Vector of persistent peaks across lambda values
#'   \item significance_intervals: Intervals showing peak significance
#'   \item lambdas: Vector of lambda values used
#'   \item time_grid: Time grid used
#'   \item mean_function: Mean function at optimal lambda
#'   \item peak_locs: Peak locations at optimal lambda
#'   \item valley_locs: Valley locations at optimal lambda
#'   \item labels: Peak labels at optimal lambda
#'   \item lambda_opt: Optimal lambda value
#'   \item num_peaks: Number of persistent peaks
#'   \item aligned_functions: Aligned functions at optimal lambda
#'   \item warping_functions: Warping functions at optimal lambda
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'my_curves' contains functional data and 't_grid' is defined
#' result <- peak_persistance_diagram(my_curves, t_grid)
#' plot(result$bc)
#' plot(result$surface)
#' }
#'
#' @importFrom graphics plot
#' @export
peak_persistance_diagram <- function(curves, t_grid,
                                     parallel = TRUE,
                                     lambda_search_start = 2,
                                     lambda_search_min_bound = 0.01,
                                     lambda_search_threshold = 1e-3,
                                     max_lambda_search_steps = 20,
                                     max_iter = 20,
                                     penalty = "roughness",
                                     n_lambda = 10,
                                     sig_threshold = 0.03,
                                     pers_threshold = 0.28,
                                     pers_method = "clustering") {

  #  Input validation:
  # # TODO:  "missing" checks for args without defaults are not necessary at all, please rm
  if (missing(curves)) {
    stop("'curves' is required")
  }
  if (!inherits(curves, "tf")) {
    stop("'curves' must be a tf object from the tidyfun package")
  }
  if (missing(t_grid)) {
    stop("'t_grid' is required")
  }
  if (!is.numeric(t_grid)) {
    stop("'t_grid' must be a numeric vector")
    # TODO: useless / not sringent enough. has to be a suitable arg vector.
    #       maybe leave check to tfb/tfd calls below
  }
  if (!is.numeric(lambda_search_start) || lambda_search_start <= 0) {
    stop("'lambda_search_start' must be a positive numeric value")
  } # use checkmate::assert_number
  if (!is.numeric(lambda_search_min_bound) || lambda_search_min_bound <= 0) {
    stop("'lambda_search_min_bound' must be a positive numeric value")
  }
  if (!is.numeric(lambda_search_threshold) || lambda_search_threshold <= 0) {
    stop("'lambda_search_threshold' must be a positive numeric value")
  }
  if (!is.numeric(max_lambda_search_steps) || max_lambda_search_steps <= 0 ||
      max_lambda_search_steps != round(max_lambda_search_steps)) {
    stop("'max_lambda_search_steps' must be a positive integer")
  }
  if (!is.numeric(max_iter) || max_iter <= 0 || max_iter != round(max_iter)) {
    stop("'max_iter' must be a positive integer")
  }
  if (!penalty %in% c("roughness", "geodesic", "norm")) {
    stop("penalty must be either 'roughness', 'geodesic', or 'norm'.")
  }
  # # TODO:  much better pattern:
  # use arg-list in function definition above, then use match.arg
  # makes docs easier to understand, tab completion works better, etc

  if (!is.numeric(n_lambda) || n_lambda <= 1 || n_lambda != round(n_lambda)) {
    stop("'n_lambda' must be an integer greater than 1")
  }
  if (!is.numeric(sig_threshold) || sig_threshold <= 0) {
    stop("'sig_threshold' must be a positive numeric value")
  }
  if (!is.numeric(pers_threshold) || pers_threshold <= 0 || pers_threshold >= 1) {
    stop("'pers_threshold' must be a positive numeric value between 0 and 1.")
  }
  if (!pers_method %in% c("clustering", "threshold")) {
    stop("'pers_method' must be either 'clustering' or 'threshold'")
  } # TODO: much better pattern: use arg-list in function definition above, then use match.arg

  # Create dataframe with function curves
  tryCatch({
    curves <- tfb(curves, basis = "spline", bs = "bs")
    # TODO: users should be able to set k (basis dimension!) for this
    # ??: shouldn't this use t_grid?
    # TODO: bad pattern : don't overwrite function inputs with new objects of same name,
    #       makes code hard to read
  }, error = function(e) {
    stop("Failed to convert curves using tfb: ", e$message)
  })

  # Find max lambda value
  # Needs tf in data representation to compute norms as this is returned by align_functions
  message("\nFinding max lambda value")
  tryCatch({
    curves_d <- tfd(curves, t_grid)
    # ??: does this always yield regular tfds (& if not do we need to guarantee this)?
  }, error = function(e) {
    stop("Failed to convert curves to tfd format: ", e$message)
  })

  lam_stop <- tryCatch({
    find_max_lambda(curves_d,
                    start_val = lambda_search_start,
                    threshold = lambda_search_threshold,
                    min_bound = lambda_search_min_bound,
                    max_iter = max_iter,
                    penalty = penalty,
                    max_search_steps = max_lambda_search_steps,
                    parallel = parallel)
  }, error = function(e) {
    stop("Error in find_max_lambda: ", e$message)
  })

  message(paste("Max lambda value: ", lam_stop))

  # Grid of lambda values
  #??: is linear lambda-grid most sensible? what about log- or sqrt-spacing?
  lambda_values <- seq(0, lam_stop, length.out = n_lambda)

  # Create column names that incorporate lambda values
  # Format: "aligned_lambda_X" where X is the lambda value
  col_names <- paste0("aligned_lambda_", formatC(lambda_values, format="f", digits=3))

  # Apply align_functions for each lambda and combine results
  max_iter_current <- max_iter
  fun_align_result_list <- lapply(lambda_values, function(lambda) {
    align_result <- tryCatch({
      align_functions(curves, lambda = lambda, parallel = parallel,
                      max_iter = max_iter_current, penalty = penalty,
                      t_grid = t_grid)
    }, error = function(e) {
      stop("Error in align_functions with lambda=", lambda, ": ", e$message)
    })

    if (align_result$converged) {
      return(list(
        functions_aligned = align_result$aligned_curves,
        warping_functions = align_result$warping_functions
      ))
    }

    if (max_iter_current == max_iter) {
      max_iter_current <<- max_iter * 2
      message("Maximum iterations doubled and trying again.")

      align_result <- tryCatch({
        align_functions(curves, lambda = lambda, parallel = parallel,
                        max_iter = max_iter_current,
                        penalty = penalty, t_grid = t_grid)
      }, error = function(e) {
        stop("Error in align_functions with doubled max_iter: ", e$message)
      })

      if (align_result$converged) {
        return(list(
          functions_aligned = align_result$aligned_curves,
          warping_functions = align_result$warping_functions
        ))
      }
    }

    message(sprintf("No convergence after %s iterations. Set max_iter higher.",
                    max_iter_current))
    return(list(
      functions_aligned = align_result$aligned_curves,
      warping_functions = align_result$warping_functions
    ))
  })

  aligned_functions_list <- lapply(fun_align_result_list, function(x) x$functions_aligned)
  warping_functions_list <- lapply(fun_align_result_list, function(x) x$warping_functions)

  # Convert to dataframe
  aligned_df <- tryCatch({
    as.data.frame(aligned_functions_list)
  }, error = function(e) {
    stop("Failed to convert aligned functions to dataframe: ", e$message)
  })
  colnames(aligned_df) <- col_names

  # Convert to dataframe
  warping_functions_df <- tryCatch({
    as.data.frame(warping_functions_list)
  }, error = function(e) {
    stop("Failed to convert warping functions to dataframe: ", e$message)
  })
  colnames(warping_functions_df) <- col_names

  # Compute shape information for each set of aligned functions
  all_peaks <- list()
  all_valleys <- list()
  mean_functions <- list()
  peak_curvatures <- list()
  peak_heights <- list()

  for(col_name in colnames(aligned_df)) {
    # Extract curves for this lambda value
    curves <- aligned_df[[col_name]]

    # Compute mean functions, find extrema and curvature and height of peaks
    results <- tryCatch({
      analyze_mean_function(curves)
    }, error = function(e) {
      stop("Error in analyze_mean_function for ", col_name, ": ", e$message)
    })

    # Store results in lists with column name as identifier
    all_peaks[[col_name]] <- results$peaks
    all_valleys[[col_name]] <- results$valleys
    mean_functions[[col_name]] <- results$function_mean
    peak_curvatures[[col_name]] <- results$peak_curvatures
    peak_heights[[col_name]] <- results$peak_heights
  }

  # Labelling peaks
  all_labels <- list()
  all_labels[[col_names[1]]] <- seq_along(all_peaks[[col_names[1]]])
  label_max <- length(all_labels[[col_names[1]]])

  # Iterate over consecutive pairs of columns, starting from second column
  for(i in 2:length(all_peaks)) {
    # Get peaks and valleys for current and previous curves
    previous_name <- col_names[[i-1]]
    curr_name <- col_names[[i]]
    peaks_1 <- all_peaks[[previous_name]]
    valleys_1 <- all_valleys[[previous_name]]
    peaks_2 <- all_peaks[[curr_name]]
    valleys_2 <- all_valleys[[curr_name]]

    # Get labels from previous curve
    labels_1 <- all_labels[[previous_name]]

    # Get new labels using peak_successor
    result <- tryCatch({
      peak_successor(
        peaks_1 = peaks_1,
        valleys_1 = valleys_1,
        peaks_2 = peaks_2,
        labels_1 = labels_1,
        label_max = label_max,
        t_grid = t_grid
      )
    }, error = function(e) {
      stop("Error in peak_successor for iteration ", i, ": ", e$message)
    })

    # Store new labels
    all_labels[[curr_name]] <- result$labels
    label_max <- result$label_max
  }

  # Calculate which peaks are significant
  significant_peaks <- lapply(col_names, function(name) {
    tryCatch({
      get_significant_peaks(peak_curvatures[[name]], all_labels[[name]], sig_threshold)
    }, error = function(e) {
      stop("Error in get_significant_peaks for ", name, ": ", e$message)
    })
  })
  names(significant_peaks) <- col_names

  # Compute peak intervals
  peak_intervals <- tryCatch({
    get_peak_intervals(significant_peaks)
  }, error = function(e) {
    stop("Error in get_peak_intervals: ", e$message)
  })

  # Calculate persistence of peaks and number of persistent peaks
  persistent_peaks <- if (pers_method == "clustering") {
    tryCatch({
      get_persistent_peaks_clustering(significant_peaks)
    }, error = function(e) {
      stop("Error in get_persistent_peaks_clustering: ", e$message)
    })
  } else {
    tryCatch({
      get_persistent_peaks_threshold(peak_intervals, pers_threshold)
    }, error = function(e) {
      stop("Error in get_persistent_peaks_threshold: ", e$message)
    })
  }

  num_peaks <- length(persistent_peaks)

  # Find optimal lambda index
  idx_opt <- tryCatch({
    find_optimal_lambda(significant_peaks, persistent_peaks)
  }, error = function(e) {
    stop("Error in find_optimal_lambda: ", e$message)
  })

  opt_lambda <- lambda_values[idx_opt]

  # Persistence diagram in barchart form
  persistence_plot_bc <- tryCatch({
    create_persistence_diagram_bc(peak_intervals)
  }, error = function(e) {
    stop("Error in create_persistence_diagram_bc: ", e$message)
  })

  # Peak Persistance surface
  persistance_plot_sf <- tryCatch({
    create_persistance_diagram_sf(t_grid, lambda_values,
                                  mean_functions, all_peaks,
                                  all_labels, significant_peaks)
  }, error = function(e) {
    stop("Error in create_persistance_diagram_sf: ", e$message)
  })

  return(list(
    bc = persistence_plot_bc,
    surface = persistance_plot_sf,
    significant_peaks = significant_peaks[[idx_opt]],
    persistent_peaks = persistent_peaks,
    significance_intervals = peak_intervals,
    lambdas = lambda_values,
    time_grid = t_grid,
    mean_function = mean_functions[[idx_opt]],
    peak_locs = all_peaks[[idx_opt]],
    valley_locs = all_valleys[[idx_opt]],
    labels = all_labels[[idx_opt]],
    lambda_opt = opt_lambda,
    num_peaks = num_peaks,
    aligned_functions = aligned_df[[col_names[[idx_opt]]]],
    warping_functions = warping_functions_df[[col_names[[idx_opt]]]]
  ))
}
