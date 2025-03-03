#' Run Single Benchmark Simulation
#'
#' @description
#' Executes a single benchmark run using the provided parameters.
#'
#' @param params A row from the benchmark grid, containing parameter settings.
#' @param function_list List of functions indexed by g_id to use for data generation.
#' @param log_file Character, path to file for logging benchmark progress.
#' @param save_plots Logical, whether to save plots generated during benchmarking.
#' @param plot_dir Character, directory to save plots if save_plots is TRUE.
#'
#' @importFrom tidyfun geom_spaghetti
#'
#' @return A list containing benchmark results, including all relevant data and results.
run_simulation <- function(params, function_list, log_file = NULL,
                           save_plots = FALSE, plot_dir = NULL,
                           seed = NULL) {
  # Initialize logging
  if (!is.null(log_file)) {
    log_msg <- paste0("Starting benchmark ID ", params$benchmark_id, " - ", format(Sys.time()))
    cat(log_msg, "\n", file = log_file, append = TRUE)
  }

  # Initialize result storage
  result <- list(
    params = params,
    runtime = list(),
    results = list()
  )

  # Start timer
  total_start_time <- Sys.time()

  # Step 1: Generate test data
  data_gen_start <- Sys.time()

  tryCatch({
    # Check if g_id is provided in params
    if (!("g_id" %in% names(params))) {
      stop("Missing g_id parameter to select a function from function_list")
    }

    # Get the generator function from the function list
    g_id <- params$g_id

    # Check valid index range:
    g_id_numeric <- as.integer(g_id)
    if (is.na(g_id_numeric) || g_id_numeric < 1 || g_id_numeric > length(function_list)) {
      stop(paste("Invalid g_id:", g_id, "- must be a number between 1 and", length(function_list)))
    }

    # Initialize parameters for data generation
    gen_params <- list()

    # Get the function:
    gen_params$g <- function_list[[g_id_numeric]]

    # Check which parameters from params should be passed to generate_functional_curves
    func_params <- names(formals(generate_functional_curves))
    for (param in names(params)) {
      if (param %in% func_params) {
        # Only include non-NA values
        if (!is.na(params[[param]])) {
          gen_params[[param]] <- params[[param]]
        }
      }
    }

    # If passed, set seed for data generation process
    if (!is.null(seed)){
      gen_params[["seed"]] <- seed
    }


    # Generate the curves
    curve_data <- do.call(generate_functional_curves, gen_params)

    # Store data generation time
    result$runtime$data_generation <- difftime(Sys.time(), data_gen_start, units = "secs")

    # Store simulated data
    result$simulated_data$curves <- curve_data$curves
    result$simulated_data$warped_grid <- curve_data$grid_warped

    # Create and store ground truth function
    base_function <- curve_data$base_function
    t_grid <- curve_data$t_grid
    ground_truth <- tfd(base_function(t_grid), t_grid)
    result$simulated_data$ground_truth <- ground_truth
    result$simulated_data$t_grid <- t_grid
    result$simulated_data$base_function <- base_function


  }, error = function(e) {
    if (!is.null(log_file)) {
      error_msg <- paste0("ERROR in data generation: ", e$message)
      cat(error_msg, "\n", file = log_file, append = TRUE)
    }
    result$errors$data_generation <- e$message
    return(result)  # Return early if data generation fails
  })

  # Step 2: Run peak persistence analysis
  if (!is.null(result$simulated_data)) {
    ppd_start <- Sys.time()

    tryCatch({
      # Extract parameters for peak_persistance_diagram
      ppd_params <- list()

      # Add curves and t_grid from curve data
      ppd_params$curves <- result$simulated_data$curves
      ppd_params$t_grid <- result$simulated_data$t_grid

      # Check which parameters from params should be passed to peak_persistance_diagram
      ppd_func_params <- names(formals(peak_persistance_diagram))
      for (param in names(params)) {
        if (param %in% ppd_func_params && param != "curves" && param != "t_grid") {
          ppd_params[[param]] <- params[[param]]
        }
      }

      # Run peak persistence analysis
      ppd_result <- do.call(peak_persistance_diagram, ppd_params)

      # Store PPD runtime
      result$runtime$peak_persistence <- difftime(Sys.time(), ppd_start, units = "secs")

      # Store the PPD results
      result$ppd_result <- ppd_result

      # Extract specific elements for easier access
      # result$results$mean_function <- ppd_result$mean_function
      # result$results$warping_functions <- ppd_result$warping_functions
      # result$results$aligned_curves <- ppd_result$aligned_functions
      # result$results$peak_locs <- ppd_result$peak_locs
      # result$results$valley_locs <- ppd_result$valley_locs
      # result$results$num_peaks <- ppd_result$num_peaks
      # result$results$persistent_peaks <- ppd_result$persistent_peaks
      # result$results$lambda_opt <- ppd_result$lambda_opt

      # Store plots
      result$plots$barchart <- ppd_result$bc
      result$plots$surface <- ppd_result$surface

    }, error = function(e) {
      if (!is.null(log_file)) {
        error_msg <- paste0("ERROR in peak persistence analysis: ", e$message)
        cat(error_msg, "\n", file = log_file, append = TRUE)
      }
      result$errors$peak_persistence <- e$message
    })
  }

  # Step 3: Run shape constrained estimation
  if (!is.null(result$ppd_result)) {
    sce_start <- Sys.time()


    tryCatch({
      # Extract parameters for shape_constrained_estimation
      sce_params <- list(
        curve_data = result$ppd_result$aligned_functions,
        peak_locs = result$ppd_result$peak_locs,
        valley_locs = result$ppd_result$valley_locs,
        significant_peaks = result$ppd_result$significant_peaks,
        peak_labels = result$ppd_result$labels,
        mean_function = result$ppd_result$mean_function,
        t_grid = result$ppd_result$time_grid
      )

      # Check which parameters from params should be passed to shape_constrained_estimation
      sce_func_params <- names(formals(shape_constrained_estimation))
      for (param in names(params)) {
        if (param %in% sce_func_params) {
          sce_params[[param]] <- params[[param]]
        }
      }
      # print(sce_params)

      # # Add this before do.call(shape_constrained_estimation, sce_params)
      # for (param_name in names(sce_params)) {
      #   message("Parameter: ", param_name, " - Class: ", paste(class(sce_params[[param_name]]), collapse=", "))
      # }

      # Run shape constrained estimation
      fn_est <- do.call(shape_constrained_estimation, sce_params)

      # Store SCE runtime
      result$runtime$shape_constrained_est <- difftime(Sys.time(), sce_start, units = "secs")

      # Store the estimated function
      result$estimated_function <- fn_est

      # Plot mean, estimated, and ground truth curve
      vis_df <- data.frame(result$ppd_result$aligned_functions)
      colnames(vis_df) <- c("curves_aligned")
      ground_truth <- result$simulated_data$ground_truth
      truth_df <- data.frame(result$simulated_data$ground_truth)
      colnames(truth_df) <- c("ground_truth")
      fn_est <- result$estimated_function
      est_df <- data.frame(fn_est)
      colnames(est_df) <- c("fn_est")
      mfn <- result$ppd_result$mean_function
      mean_df <- data.frame(mfn)
      colnames(mean_df) <- c("mfn")

      plot_fn <- ggplot() +
        geom_spaghetti(data = vis_df, aes(y = curves_aligned), alpha = 0.3) +
        geom_spaghetti(data = est_df, aes(y = fn_est), color = "orange", linewidth = 1.5) +
        geom_spaghetti(data = mean_df, aes(y = mfn), color = "green", linewidth = 1.5) +
        geom_spaghetti(data = truth_df, aes(y = ground_truth), color = "red", linewidth = 1.5) +
        labs(caption = "red: true, green: mean, orange: shape-constrained")

      result$plots$estimates <- plot_fn

      #TODO: maybe plot warping functions next to each other as well

    }, error = function(e) {
      if (!is.null(log_file)) {
        error_msg <- paste0("ERROR in shape constrained estimation: ", e$message)
        cat(error_msg, "\n", file = log_file, append = TRUE)
      }
      result$errors$shape_constrained <- e$message
    })
  }

  # Step 4: Align original curves to estimated function
  if (!is.null(result$estimated_function) && !is.null(result$simulated_data$curves)) {
    align_start <- Sys.time()

    tryCatch({
      # Extract parameters
      original_curves <- result$simulated_data$curves
      estimated_function <- result$estimated_function
      t_grid <- result$simulated_data$t_grid

      # Default lambda value
      lambda_align <- 0.1
      if ("lambda_align" %in% names(params)) {
        lambda_align <- params$lambda_align
      }

      # Align original curves to estimated function
      alignment_result <- align_to_mean(
        function_curves = original_curves,
        target_function = estimated_function,
        t_grid = t_grid,
        lambda = lambda_align
      )

      # Store alignment runtime
      result$runtime$alignment <- difftime(Sys.time(), align_start, units = "secs")

      # Store alignment results
      result$alignment_to_est <- alignment_result$aligned_functions
      result$al_to_est_warping_functions <- alignment_result$warping_functions

    }, error = function(e) {
      if (!is.null(log_file)) {
        error_msg <- paste0("ERROR in alignment to estimated function: ", e$message)
        cat(error_msg, "\n", file = log_file, append = TRUE)
      }
      result$errors$alignment <- e$message
    })
  }

  # Save plots if requested
  if (save_plots && !is.null(plot_dir) && !is.null(result$plots)) {
    tryCatch({
      plot_dir_path <- file.path(plot_dir, paste0("benchmark_", params$benchmark_id))
      dir.create(plot_dir_path, recursive = TRUE, showWarnings = FALSE)

      # Save PPD plots if available
      if (!is.null(result$plots$barchart)) {
        ggsave(file.path(plot_dir_path, "barchart.png"), result$plots$barchart)
      }
      if (!is.null(result$plots$surface)) {
        ggsave(file.path(plot_dir_path, "surface.png"), result$plots$surface)
      }

      # Save original curves plot if available
      if (!is.null(result$simulated_data)) {
        original_plot <- plot_simulated_curves(
          result$simulated_data$curves,
          result$simulated_data$t_grid,
          result$simulated_data$base_function
        )
        ggsave(file.path(plot_dir_path, "original_curves.png"), original_plot)
      }

      if (!is.null(result$plots$estimates)) {
        ggsave(file.path(plot_dir_path, "estimates.png"), result$plots$estimates)
      }

    }, error = function(e) {
      if (!is.null(log_file)) {
        error_msg <- paste0("ERROR in saving plots: ", e$message)
        cat(error_msg, "\n", file = log_file, append = TRUE)
      }
      result$errors$plots <- e$message
    })
  }

  # Calculate total runtime
  result$runtime$total <- difftime(Sys.time(), total_start_time, units = "secs")

  # Log completion
  if (!is.null(log_file)) {
    completion_msg <- paste0(
      "Completed benchmark ID ", params$benchmark_id,
      " in ", round(result$runtime$total, 2), " seconds - ",
      format(Sys.time())
    )
    cat(completion_msg, "\n", file = log_file, append = TRUE)
  }

  return(result)
}
