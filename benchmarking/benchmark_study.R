#' Run Benchmark Study
#'
#' @description
#' Executes a complete benchmark study by running multiple benchmarks with different
#' parameter combinations and calculating metrics. Saves results after each run to
#' prevent data loss in case of crashes.
#'
#' @param param_grid Data frame containing parameter combinations to benchmark.
#' @param function_list List of functions to use for data generation.
#' @param output_dir Character, directory for saving results.
#' @param save_plots Logical, whether to save plots generated during benchmarking.
#' @param parallel_processing Logical, whether to use parallel processing.
#' @param n_cores Integer, number of cores to use if parallel_processing is TRUE.
#'
#' @return A list containing all benchmark results and metrics.
#'
#' @examples
#' \dontrun{
#' # Create parameter grid
#' grid <- create_benchmark_grid(
#'   g_id = c(1, 2),
#'   n_lambda = c(5, 10),
#'   lambda_search_min_bound = c(0.1, 0.5),
#'   warping = list(
#'     simple = list(),
#'     flexible = list(
#'       warping_gamma = c(1, 5),
#'       warping_points = c(3, 5)
#'     )
#'   )
#' )
#'
#' # Create function list
#' function_list <- list(func_1, func_2)
#'
#' # Run benchmark study
#' study_results <- run_benchmark_study(
#'   param_grid = grid,
#'   function_list = function_list,
#'   output_dir = "./benchmark_results",
#'   save_plots = TRUE
#' )
#' }
#'
#' @export
run_benchmark_study <- function(param_grid, function_list,
                                output_dir = "./benchmarking/results/test",
                                save_plots = TRUE, parallel_processing = FALSE,
                                n_cores = 2, seed = 123) {
  # Input validation
  if (!is.data.frame(param_grid)) {
    stop("param_grid must be a data frame")
  }

  if (!is.list(function_list)) {
    stop("function_list must be a list of functions")
  }

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Initialize log file
  log_file <- file.path(output_dir, "benchmark_study_log.txt")
  cat("Starting benchmark study: ", format(Sys.time()), "\n", file = log_file)
  cat("Number of parameter combinations: ", nrow(param_grid), "\n", file = log_file, append = TRUE)

  # Create metrics dataframe for incremental saving
  metrics_df_file <- file.path(output_dir, "metrics_dataframe.rds")
  results_dir <- file.path(output_dir, "individual_results")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # Check if a previous run exists and load the metrics dataframe
  if (file.exists(metrics_df_file)) {
    metrics_df <- readRDS(metrics_df_file)
    completed_ids <- metrics_df$benchmark_id
    cat("Found existing metrics dataframe with", length(completed_ids), "completed benchmarks\n",
        file = log_file, append = TRUE)
  } else {
    metrics_df <- data.frame()
    completed_ids <- c()
  }

  # # Check for parallel package if using parallel processing
  # if (parallel_processing) {
  #   if (!requireNamespace("parallel", quietly = TRUE)) {
  #     warning("'parallel' package is needed for parallel processing. Falling back to sequential execution.")
  #     parallel_processing <- FALSE
  #   }
  # }

  # Create plot directory if needed
  plot_dir <- NULL
  if (save_plots) {
    plot_dir <- file.path(output_dir, "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Function to convert a metrics list to a data frame row
  metrics_to_df_row <- function(metrics, benchmark_id) {
    # Start with the benchmark ID
    row_data <- list(benchmark_id = benchmark_id)

    # Add all parameters
    params <- metrics$params
    for (param_name in names(params)) {
      row_data[[param_name]] <- params[[param_name]]
    }

    # Add key metrics
    add_metric <- function(metric_name, prefix = "") {
      parts <- strsplit(metric_name, "$", fixed = TRUE)[[1]]
      value <- metrics

      for (part in parts) {
        if (is.null(value) || !part %in% names(value)) {
          return(NA)
        }
        value <- value[[part]]
      }

      col_name <- paste0(prefix, gsub("$", "_", metric_name, fixed = TRUE))
      row_data[[col_name]] <- value
    }

    # Add runtime metrics
    if (!is.null(metrics$runtime)) {
      for (rt_name in names(metrics$runtime)) {
        row_data[[paste0("runtime_", rt_name)]] <- as.numeric(metrics$runtime[[rt_name]])
      }
    }

    # Add distances
    add_metric("mean_ground_truth_distance")
    add_metric("est_ground_truth_distance")

    # Add warping distances
    if (!is.null(metrics$warping_distances)) {
      for (wd_name in names(metrics$warping_distances)) {
        row_data[[paste0("warping_distances_", wd_name)]] <- metrics$warping_distances[[wd_name]]
      }
    }

    # Add peak detection metrics
    if (!is.null(metrics$peak_detection_metrics)) {
      for (pd_name in names(metrics$peak_detection_metrics)) {
        row_data[[paste0("peak_detection_", pd_name)]] <- metrics$peak_detection_metrics[[pd_name]]
      }
    }

    # Add alignment scores
    if (!is.null(metrics$alignment_scores)) {
      for (as_name in names(metrics$alignment_scores)) {
        row_data[[paste0("alignment_", as_name)]] <- metrics$alignment_scores[[as_name]]
      }
    }

    # Convert to data frame row
    return(data.frame(row_data, stringsAsFactors = FALSE))
  }

  # Function to run a single benchmark, calculate metrics, and save results
  run_single_benchmark <- function(i) {
    benchmark_id <- param_grid$benchmark_id[i]

    # Check if this benchmark has already been completed
    if (benchmark_id %in% completed_ids) {
      cat("Skipping benchmark", benchmark_id, "(already completed)\n")
      return(NULL)
    }

    cat("Running benchmark", i, "of", nrow(param_grid), "(ID:", benchmark_id, ")\n")

    # Get parameters for this run
    params <- param_grid[i, ]

    # Run the benchmark
    tryCatch({
      benchmark_result <- run_simulation(
        params = params,
        function_list = function_list,
        log_file = log_file,
        save_plots = save_plots,
        plot_dir = plot_dir,
        seed = seed
      )

      # Calculate metrics
      metrics <- calculate_benchmark_metrics(benchmark_result)

      # Save individual result
      saveRDS(benchmark_result, file.path(results_dir, paste0("result_", benchmark_id, ".rds")))
      saveRDS(metrics, file.path(results_dir, paste0("metrics_", benchmark_id, ".rds")))

      # Convert metrics to data frame row
      metrics_row <- metrics_to_df_row(metrics, benchmark_id)

      # Return the benchmark result and metrics
      return(list(
        result = benchmark_result,
        metrics = metrics,
        metrics_row = metrics_row
      ))
    }, error = function(e) {
      # Log error and continue with next benchmark
      cat("ERROR in benchmark", benchmark_id, ":", e$message, "\n", file = log_file, append = TRUE)
      return(NULL)
    })
  }

  # Run the benchmarks (sequentially to ensure reliable incremental saving)
  all_results <- list()

  for (i in 1:nrow(param_grid)) {
    result <- run_single_benchmark(i)

    if (!is.null(result)) {
      # Add to results list
      all_results[[length(all_results) + 1]] <- result

      # Update metrics dataframe
      metrics_df <- rbind(metrics_df, result$metrics_row)

      # Save updated metrics dataframe after each successful run
      saveRDS(metrics_df, metrics_df_file)

      # Update log
      cat("Completed benchmark", param_grid$benchmark_id[i], "\n", file = log_file, append = TRUE)
    }
  }

  # Extract final results
  benchmark_results <- lapply(all_results, function(r) r$result)
  all_metrics <- lapply(all_results, function(r) r$metrics)

  # Save final consolidated results
  saveRDS(benchmark_results, file.path(output_dir, "all_benchmark_results.rds"))
  saveRDS(all_metrics, file.path(output_dir, "all_metrics.rds"))

  # Save metrics dataframe as CSV for easier analysis
  if (nrow(metrics_df) > 0) {
    utils::write.csv(metrics_df, file.path(output_dir, "metrics_dataframe.csv"), row.names = FALSE)
  }

  # Update log file
  cat("Benchmark study completed: ", format(Sys.time()), "\n", file = log_file, append = TRUE)
  cat("Results saved to: ", output_dir, "\n", file = log_file, append = TRUE)
  cat("Total benchmarks completed: ", nrow(metrics_df), " of ", nrow(param_grid), "\n",
      file = log_file, append = TRUE)

  # Return results
  return(list(
    benchmark_results = benchmark_results,
    metrics = all_metrics,
    metrics_df = metrics_df,
    param_grid = param_grid
  ))
}
