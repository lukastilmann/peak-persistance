library(fdasrvf)

# creates barchart and surface diagrams tracking peaks of aligned functions
# for different lambda values
peak_persistance_diagram <- function(curves, t_grid, lambda_search_start = 2,
                                     lambda_search_min_bound = 0.01,
                                     lambda_search_threshold = 1e-3,
                                     n_lambda = 10,
                                     sig_threshold = 0.03, pers_threshold = 0.28,
                                     align_fun = "tw"){
  # Create dataframe with function curves
  curves <- tfb(curves, basis = "spline", bs = "bs")
  curves_df <- data.frame(curves)
  colnames(curves_df) <- c("curves")

  # Find max lambda value
  # Needs tf in data representation to compute norms as this is returned by align_functions
  print("Finding max lambda value")
  curves_d <- tfd(curves, t_grid)
  curves_d_df <- data.frame(curves_d)
  colnames(curves_d_df) <- c("curves")
  lam_stop <- find_max_lambda(curves_d_df, t_grid,
                              start_val = lambda_search_start,
                              threshold = lambda_search_threshold,
                              max_iter = 20,
                              min_bound = lambda_search_min_bound,
                              parallel = TRUE)
  print(paste("Max lambda value: ", lam_stop))

  # Grid of lambda values
  lambda_values <- seq(0, lam_stop, length.out = n_lambda)

  # Create column names that incorporate lambda values
  # Format: "aligned_lambda_X" where X is the lambda value
  col_names <- paste0("aligned_lambda_", formatC(lambda_values, format="f", digits=3))

  # Apply align_functions for each lambda and combine results
  aligned_list <- lapply(lambda_values, function(lambda) {
    print(paste("Alignment for:", lambda))
    aligned <- align_functions(curves_df, lambda = lambda, parallel = TRUE,
                               max_iter = 20, t_grid = t_grid,
                               func = align_fun)
    return(aligned)
  })

  # Convert to dataframe
  aligned_df<- as.data.frame(aligned_list)
  colnames(aligned_df) <- col_names

  # Add lambda values as attribute for easy access later
  #attr(aligned_df, "lambda_values") <- lambda_values

  # Compute information about mean function shape
  all_peaks <- list()
  all_valleys <- list()
  mean_functions <- list()
  peak_curvatures <- list()
  peak_heights <- list()

  # For each column
  for(col_name in colnames(aligned_df)) {
    # Extract curves for this lambda value
    curves <- aligned_df[[col_name]]

    # Compute mean functions, find extrema and curvature and height of peaks
    results <- analyze_mean_function(curves)

    # Store results in lists with column name as identifier
    all_peaks[[col_name]] <- results$peaks
    all_valleys[[col_name]] <- results$valleys
    mean_functions[[col_name]] <- results$function_mean
    peak_curvatures[[col_name]] <- results$peak_curvatures
    peak_heights[[col_name]] <- results$peak_heights
  }

  # Labelling peaks
  all_labels <- list()
  # Labels for first mean function
  all_labels[[col_names[1]]] <- seq_along(all_peaks[[col_names[1]]])
  max_label = length(all_labels[[col_names[1]]])
  print(max_label)

  # Iterate over consecutive pairs of columns, starting from second column
  for(i in 2:length(all_peaks)) {
    # Get peaks and valleys for current and previous curves
    previous_name <- col_names[[i-1]]
    curr_name <- col_names[[i]]

    peaks1 <- all_peaks[[previous_name]]
    valleys1 <- all_valleys[[previous_name]]
    peaks2 <- all_peaks[[curr_name]]
    valleys2 <- all_valleys[[curr_name]]

    # Get labels from previous curve
    labels1 <- all_labels[[previous_name]]

    # Get new labels using peak_successor
    result <- peak_successor(
      peaks1 = peaks1,
      valleys1 = valleys1,
      peaks2 = peaks2,
      valleys2 = valleys2,
      labels1 = labels1,
      label_max = max_label,
      t_grid = t_grid
    )

    # Store new labels
    all_labels[[curr_name]] <- result$labels
    max_label = result$label_max
  }

  # Calculate which peaks are significant
  significant_peaks <- lapply(col_names, function(name){
    get_significant_peaks(peak_curvatures[[name]], all_labels[[name]], sig_threshold)
  })

  names(significant_peaks) <- col_names


  # Calculate persistance of peaks
  res <- find_persistent_peaks(significant_peaks, pers_threshold)
  persistent_peaks <- res$persistent_labels
  num_peaks <- length(res$persistent_labels)

  # Find optimal lambda index
  idx_opt <- find_optimal_lambda(significant_peaks,
                                 persistent_peaks)
  opt_lambda <- lambda_values[idx_opt]

  # Persistance diagram in barchart form
  persistence_plot_bc <- create_persistence_diagram_bc(res$intervals)

  # Peak Persistance surface
  persistance_plot_sf <- create_persistance_diagram_sf(t_grid, lambda_values,
                                                       mean_functions, all_peaks,
                                                       all_labels, significant_peaks)

  return(list(
    bc = persistence_plot_bc,
    surface = persistance_plot_sf,
    significant_peaks = significant_peaks,
    persistent_peaks = res$persistent_labels,
    intervals = res$intervals,
    lam = lambda_values,
    grid = t_grid,
    mfn = mean_functions,
    peak_locs = all_peaks,
    valley_locs = all_valleys,
    labels = all_labels,
    opt_lam_idx = idx_opt,
    opt_lam = opt_lambda,
    num_peaks = num_peaks,
    aligned_functions = aligned_df
  ))
}
