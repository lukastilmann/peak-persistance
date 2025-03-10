#' Run Benchmark Study
#'
#' @description
#' Executes a complete benchmark study by running multiple benchmarks with different
#' parameter combinations and calculating metrics. Saves results after each run to
#' prevent data loss in case of crashes, while minimizing storage requirements.
#'
#' @param param_grid Data frame containing parameter combinations to benchmark.
#' @param function_list List of functions to use for data generation.
#' @param output_dir Character, directory for saving results.
#' @param save_plots Logical, whether to save plots generated during benchmarking.
#' @param parallel_processing Logical, whether to use parallel processing.
#' @param n_cores Integer, number of cores to use if parallel_processing is TRUE.
#' @param seed Integer or NULL, base seed for random number generation.
#' @param runs_per_config Integer, number of repetitions to run for each parameter combination.
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
#'   save_plots = TRUE,
#'   runs_per_config = 5
#' )
#' }
#'
#' @export
run_benchmark_study <- function(param_grid, function_list,
                                output_dir = "./benchmarking/results/test",
                                save_plots = TRUE, seed = NULL,
                                runs_per_config = 5) {
  # Input validation
  if (!is.data.frame(param_grid)) {
    stop("param_grid must be a data frame")
  }
  if (!is.list(function_list)) {
    stop("function_list must be a list of functions")
  }
  if (!is.numeric(runs_per_config) || runs_per_config < 1) {
    stop("runs_per_config must be a positive integer")
  }

  # Set base random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Initialize log file
  log_file <- file.path(output_dir, "benchmark_study_log.txt")
  cat("Starting benchmark study: ", format(Sys.time()), "\n", file = log_file)
  cat("Number of parameter combinations: ", nrow(param_grid), "\n", file = log_file, append = TRUE)
  cat("Runs per configuration: ", runs_per_config, "\n", file = log_file, append = TRUE)
  cat("Total benchmark runs: ", nrow(param_grid) * runs_per_config, "\n", file = log_file, append = TRUE)

  # Expand parameter grid to include multiple runs per configuration
  expanded_grid <- expand_runs_grid(param_grid, runs_per_config, seed)
  cat("Created expanded parameter grid with ", nrow(expanded_grid), " total runs\n",
      file = log_file, append = TRUE)

  # File paths for incremental saving
  metrics_df_file <- file.path(output_dir, "metrics_dataframe.rds")
  all_results_file <- file.path(output_dir, "all_benchmark_results.rds")
  all_metrics_file <- file.path(output_dir, "all_metrics.rds")

  # Check if a previous run exists and load the metrics dataframe and previous results
  if (file.exists(metrics_df_file)) {
    metrics_df <- readRDS(metrics_df_file)
    completed_ids <- metrics_df$benchmark_id
    cat("Found existing metrics dataframe with", length(completed_ids), "completed benchmarks\n",
        file = log_file, append = TRUE)

    # Load existing results if available
    if (file.exists(all_results_file)) {
      benchmark_results <- readRDS(all_results_file)
      cat("Loaded existing benchmark results\n", file = log_file, append = TRUE)
    } else {
      benchmark_results <- list()
    }

    # Load existing metrics if available
    if (file.exists(all_metrics_file)) {
      all_metrics <- readRDS(all_metrics_file)
      cat("Loaded existing metrics\n", file = log_file, append = TRUE)
    } else {
      all_metrics <- list()
    }
  } else {
    metrics_df <- data.frame()
    completed_ids <- c()
    benchmark_results <- list()
    all_metrics <- list()
  }

  # Create plot directory if needed
  plot_dir <- NULL
  if (save_plots) {
    plot_dir <- file.path(output_dir, "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Function to convert a metrics list to a data frame row
  metrics_to_df_row <- function(metrics, benchmark_id, run_seed, run_index) {
    # Start with the benchmark ID, run index, and seed
    row_data <- list(
      benchmark_id = benchmark_id,
      run_index = run_index,
      run_seed = run_seed
    )

    # Add all parameters
    params <- metrics$params
    for (param_name in names(params)) {
      row_data[[param_name]] <- params[[param_name]]
    }

    # Add runtime metrics
    if (!is.null(metrics$runtime)) {
      for (rt_name in names(metrics$runtime)) {
        row_data[[paste0("runtime_", rt_name)]] <- as.numeric(metrics$runtime[[rt_name]])
      }
    }

    # Add number of peaks from peak-persistence-diagram
    row_data[["num_peaks"]] <- metrics$num_peaks

    # Add distances
    row_data[["mean_ground_truth_distance"]] <- metrics$mean_ground_truth_distance
    row_data[["est_ground_truth_distance"]] <- metrics$est_ground_truth_distance

    # Add warping distances
    if (!is.null(metrics$warping_distances)) {
      for (wd_name in names(metrics$warping_distances)) {
        row_data[[paste0("warping_distances_", wd_name)]] <- metrics$warping_distances[[wd_name]]
      }
    }

    # Convert to data frame row
    return(data.frame(row_data, stringsAsFactors = FALSE))
  }

  # Function to run a single benchmark, calculate metrics, and save results
  run_single_benchmark <- function(i) {
    benchmark_id <- expanded_grid$benchmark_id[i]
    run_seed <- expanded_grid$run_seed[i]
    run_index <- expanded_grid$run_index[i]
    config_id <- expanded_grid$config_id[i]

    # Check if this benchmark has already been completed
    if (benchmark_id %in% completed_ids) {
      cat("Skipping benchmark", benchmark_id, "(already completed)\n")
      return(NULL)
    }

    cat("Running benchmark", i, "of", nrow(expanded_grid),
        "(ID:", benchmark_id, ", Config:", config_id,
        ", Run:", run_index, ", Seed:", run_seed, ")\n")

    # Get parameters for this run
    params <- expanded_grid[i, ]

    # Run the benchmark
    tryCatch({
      benchmark_result <- run_simulation(
        params = params,
        function_list = function_list,
        log_file = log_file,
        save_plots = save_plots,
        plot_dir = plot_dir,
        seed = run_seed  # Use the unique seed for this run
      )

      # Store the run seed and run index in the result
      benchmark_result$run_seed <- run_seed
      benchmark_result$run_index <- run_index
      benchmark_result$config_id <- config_id
      benchmark_result$benchmark_id <- benchmark_id

      # Calculate metrics
      metrics <- calculate_benchmark_metrics(benchmark_result)

      # Convert metrics to data frame row
      metrics_row <- metrics_to_df_row(metrics, benchmark_id, run_seed, run_index)

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
  for (i in 1:nrow(expanded_grid)) {
    result <- run_single_benchmark(i)
    if (!is.null(result)) {
      # Add to results and metrics collections
      benchmark_results[[length(benchmark_results) + 1]] <- result$result
      all_metrics[[length(all_metrics) + 1]] <- result$metrics

      # Update metrics dataframe
      metrics_df <- rbind(metrics_df, result$metrics_row)

      # Save updated files after each successful run
      saveRDS(metrics_df, metrics_df_file)
      saveRDS(benchmark_results, all_results_file)
      saveRDS(all_metrics, all_metrics_file)

      # Update log
      cat("Completed benchmark", expanded_grid$benchmark_id[i],
          ", updated checkpoint files\n", file = log_file, append = TRUE)
    }
  }

  # Save expanded grid for reference
  saveRDS(expanded_grid, file.path(output_dir, "expanded_param_grid.rds"))

  # Save metrics dataframe as CSV for easier analysis
  if (nrow(metrics_df) > 0) {
    utils::write.csv(metrics_df, file.path(output_dir, "metrics_dataframe.csv"), row.names = FALSE)
  }

  # Update log file
  cat("Benchmark study completed: ", format(Sys.time()), "\n", file = log_file, append = TRUE)
  cat("Results saved to: ", output_dir, "\n", file = log_file, append = TRUE)
  cat("Total benchmarks completed: ", nrow(metrics_df), " of ", nrow(expanded_grid), "\n",
      file = log_file, append = TRUE)

  # Return results
  return(list(
    benchmark_results = benchmark_results,
    metrics = all_metrics,
    metrics_df = metrics_df,
    param_grid = param_grid,
    expanded_grid = expanded_grid
  ))
}

#' Expand Parameter Grid for Multiple Runs
#'
#' @description
#' Takes a parameter grid and expands it to include multiple runs per parameter
#' combination, each with a unique benchmark ID and random seed.
#'
#' @param param_grid Data frame containing parameter combinations to benchmark.
#' @param runs_per_config Integer, number of repetitions to run for each parameter combination.
#' @param base_seed Integer or NULL, base seed for random number generation.
#'
#' @return A data frame with expanded parameter combinations, including run indices and seeds.
#'
#' @keywords internal
expand_runs_grid <- function(param_grid, runs_per_config, base_seed = NULL) {
  # If the grid doesn't have a config_id column, add one
  if (!"config_id" %in% names(param_grid)) {
    param_grid$config_id <- 1:nrow(param_grid)
  }

  # Generate deterministic seeds based on base_seed if provided
  if (!is.null(base_seed)) {
    set.seed(base_seed)
  }

  # Create an expanded grid with runs_per_config rows for each original row
  expanded_list <- list()
  for (i in 1:nrow(param_grid)) {
    row_data <- param_grid[i, , drop = FALSE]
    orig_config_id <- row_data$config_id
    for (j in 1:runs_per_config) {
      # Create a copy of the row
      new_row <- row_data
      # Generate a unique random seed for this run
      run_seed <- sample.int(1e6, 1)
      # Set run index
      new_row$run_index <- j
      # Set run seed
      new_row$run_seed <- run_seed
      # Create a unique benchmark_id by combining config_id and run_index
      new_row$benchmark_id <- paste0(orig_config_id, "_", j)
      # Add to the list
      expanded_list[[length(expanded_list) + 1]] <- new_row
    }
  }

  # Combine all rows into a single data frame
  expanded_grid <- do.call(rbind, expanded_list)
  rownames(expanded_grid) <- NULL
  return(expanded_grid)
}
