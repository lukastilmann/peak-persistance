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
#' @importFrom parallel detectCores makeCluster stopCluster clusterExport clusterEvalQ
#' @importFrom foreach foreach '%dopar%'
#' @importFrom doParallel registerDoParallel
#' @export
peak_persistance_diagram <- function(curves, t_grid,
                                     parallel = FALSE,
                                     parallel_max_lambda = FALSE,
                                     lambda_search_start = 2,
                                     lambda_search_min_bound = 0.01,
                                     lambda_search_threshold = 1e-3,
                                     max_lambda_search_steps = 20,
                                     curvature_percentile = NULL,
                                     basis_dim = 20,
                                     max_iter = 20,
                                     penalty = c("roughness", "geodesic"),
                                     n_lambda = 10,
                                     lambda_grid_spacing = "sqrt",
                                     sig_threshold = 0.03,
                                     pers_threshold = 0.28,
                                     pers_method = c("clustering", "threshold")) {

  #  Input validation:
  checkmate::assert_class(curves, "tf")
  checkmate::assert_number(lambda_search_start, lower = 0, finite = TRUE)
  checkmate::assert_number(lambda_search_min_bound, lower = 0, finite = TRUE)
  checkmate::assert_number(lambda_search_threshold, lower = 0, finite = TRUE)
  checkmate::assert_count(max_lambda_search_steps, positive = TRUE)
  checkmate::assert_count(max_iter, positive = TRUE)
  checkmate::assert_int(n_lambda, lower = 2)
  checkmate::assert_number(sig_threshold, lower = 0, finite = TRUE)
  checkmate::assert_number(pers_threshold, lower = 0, upper = 1, finite = TRUE)

  penalty <- match.arg(penalty)
  pers_method <- match.arg(pers_method)

  # Create dataframe with function curves
  tryCatch({
    curves_b <- tfb(curves, arg = t_grid, basis = "spline",
                    bs = "bs", k = basis_dim)
  }, error = function(e) {
    stop("Failed to convert curves using tfb: ", e$message)
  })

  # Find significance threshold if curvature_percentile passed
  if (!is.null(curvature_percentile)){
    sig_threshold <- find_tau(curves_b, t_grid, curvature_percentile)
  }

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
                    parallel = parallel_max_lambda)
  }, error = function(e) {
    stop("Error in find_max_lambda: ", e$message)
  })

  message(paste("Max lambda value: ", lam_stop))

  # Grid of lambda values
  lambda_values <- create_lambda_grid(lam_stop, n_lambda, lambda_grid_spacing)

  # Create column names that incorporate lambda values
  # Format: "aligned_lambda_X" where X is the lambda value
  col_names <- paste0("aligned_lambda_", format_lambda(lambda_values))

  # Set up parallel processing if passed as argument
  if (parallel) {
    # Check if parallel package is available
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' is not available. Using sequential processing instead.")
      parallel <- FALSE
    } else {
      # Determine number of cores to use
      n_cores <- parallel::detectCores() - 1
      if (n_cores < 1) n_cores <- 1
      if (n_cores > length(lambda_values)) n_cores <- length(lambda_values)

      # Create cluster
      cl <- parallel::makeCluster(n_cores)

      # Export necessary objects and functions to the cluster
      parallel::clusterExport(cl,
                              varlist = c("curves_b", "t_grid", "max_iter", "penalty"),
                              envir = environment())

      # Load required packages on all cluster nodes
      parallel::clusterEvalQ(cl, {
        library(tidyfun)
      })

      # Register the parallel backend with foreach
      if (requireNamespace("doParallel", quietly = TRUE)) {
        doParallel::registerDoParallel(cl)
      } else {
        warning("Package 'doParallel' is not available. Using sequential processing instead.")
        parallel <- FALSE
        parallel::stopCluster(cl)
      }
    }
  }

  # Process each lambda in parallel or sequentially
  if (parallel && requireNamespace("foreach", quietly = TRUE) &&
      requireNamespace("doParallel", quietly = TRUE)) {

    # Use foreach to parallelize across lambda values
    fun_align_result_list <- foreach::foreach(
      lambda = lambda_values,
      .export = c("align_functions")
    ) %dopar% {
      max_iter_current <- max_iter

      # First attempt with regular max_iter
      align_result <- tryCatch({
        align_functions(curves_b, lambda = lambda, parallel = FALSE,
                        max_iter = max_iter_current, penalty = penalty,
                        t_grid = t_grid)
      }, error = function(e) {
        message("Error in align_functions with lambda=", lambda, ": ", e$message)
        return(NULL)
      })

      # If not converged, try with doubled max_iter
      if (!is.null(align_result) && !align_result$converged) {
        max_iter_current <- max_iter * 2
        message("Maximum iterations doubled and trying again for lambda=", lambda)

        align_result <- tryCatch({
          align_functions(curves_b, lambda = lambda, parallel = FALSE,
                          max_iter = max_iter_current,
                          penalty = penalty, t_grid = t_grid)
        }, error = function(e) {
          message("Error in align_functions with doubled max_iter for lambda=", lambda, ": ", e$message)
          return(NULL)
        })

        if (!is.null(align_result) && !align_result$converged) {
          message(sprintf("No convergence after %s iterations for lambda=%s. Consider setting max_iter higher.",
                          max_iter_current, lambda))
        }
      }

      if (is.null(align_result)) {
        # Return empty result in case of error
        return(list(
          functions_aligned = NULL,
          warping_functions = NULL
        ))
      }

      return(list(
        functions_aligned = align_result$aligned_curves,
        warping_functions = align_result$warping_functions
      ))
    }

    # Clean up parallel resources
    parallel::stopCluster(cl)

  } else {
    # Sequential processing

    fun_align_result_list <- lapply(lambda_values, function(lambda) {
      max_iter_current <- max_iter

      align_result <- tryCatch({
        align_functions(curves_b, lambda = lambda, parallel = FALSE,
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
        max_iter_current <- max_iter * 2
        message("Maximum iterations doubled and trying again for lambda=", lambda)

        align_result <- tryCatch({
          align_functions(curves_b, lambda = lambda, parallel = FALSE,
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

      message(sprintf("No convergence after %s iterations for lambda=%s. Set max_iter higher.",
                      max_iter_current, lambda))

      return(list(
        functions_aligned = align_result$aligned_curves,
        warping_functions = align_result$warping_functions
      ))
    })
  }

  # Check for NULL results
  has_nulls <- any(sapply(fun_align_result_list, function(x) is.null(x$functions_aligned)))
  if (has_nulls) {
    warning("Some lambda values failed to produce valid results. Check the messages for details.")
  }

  # Extract aligned functions and warping functions
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

  # Parallelize the shape information analysis
  if (parallel && requireNamespace("foreach", quietly = TRUE) &&
      requireNamespace("doParallel", quietly = TRUE)) {

    # Set up parallel cluster
    n_cores <- parallel::detectCores() - 1
    if (n_cores < 1) n_cores <- 1
    if (n_cores > length(col_names)) n_cores <- length(col_names)
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    # Export the aligned dataframe to all workers
    parallel::clusterExport(cl, c("aligned_df"), envir = environment())

    # Parallel processing of shape information
    results_list <- foreach::foreach(
      col_name = col_names,
      .export = c("analyze_mean_function")
    ) %dopar% {
      # Extract curves for this lambda value
      curves_al <- aligned_df[[col_name]]

      # Compute mean functions, find extrema and curvature and height of peaks
      tryCatch({
        results <- analyze_mean_function(curves_al)
        list(
          col_name = col_name,
          peaks = results$peaks,
          valleys = results$valleys,
          function_mean = results$function_mean,
          peak_curvatures = results$peak_curvatures,
          peak_heights = results$peak_heights,
          success = TRUE
        )
      }, error = function(e) {
        list(
          col_name = col_name,
          message = paste("Error in analyze_mean_function for", col_name, ":", e$message),
          success = FALSE
        )
      })
    }

    # Clean up parallel resources
    parallel::stopCluster(cl)

    # Process results and handle any errors
    for (result in results_list) {
      if (result$success) {
        col_name <- result$col_name
        all_peaks[[col_name]] <- result$peaks
        all_valleys[[col_name]] <- result$valleys
        mean_functions[[col_name]] <- result$function_mean
        peak_curvatures[[col_name]] <- result$peak_curvatures
        peak_heights[[col_name]] <- result$peak_heights
      } else {
        stop(result$message)
      }
    }
  } else {
    # Sequential processing as fallback
    for(col_name in colnames(aligned_df)) {
      # Extract curves for this lambda value
      curves_al <- aligned_df[[col_name]]

      # Compute mean functions, find extrema and curvature and height of peaks
      results <- tryCatch({
        analyze_mean_function(curves_al)
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
    sig_threshold = sig_threshold,
    aligned_functions = aligned_df[[col_names[[idx_opt]]]],
    warping_functions = warping_functions_df[[col_names[[idx_opt]]]]
  ))
}
