library(fdasrvf)
library(tidyfun)
library(tidyverse)
library(ggplot2)


tf_to_fdasrvf <- function(tf_data, tf_col, grid = NULL) {
  # Convert tf object to wide format using tf_spread with specified grid
  wide_data <- tf_spread(tf_data, !!enquo(tf_col), arg = grid)

  # Get the column names that contain the function evaluations
  value_cols <- grep(paste0(quo_name(enquo(tf_col)), "_"), names(wide_data), value = TRUE)

  # Extract just the function values into a matrix
  func_matrix <- as.matrix(wide_data[, value_cols])

  # Get dimensions
  N <- nrow(func_matrix)  # number of curves
  T <- ncol(func_matrix)  # number of time points
  n <- 2  # dimension (R²)

  # Create the 3D array with proper dimensions: n × T × N
  beta <- array(0, dim = c(n, T, N))

  # Fill the array:
  # First coordinate (x) is the grid points
  # Second coordinate (y) is the function values
  beta[1, , ] <- matrix(grid, nrow = 1, ncol = T) %>%
    replicate(n = N, simplify = "array")
  beta[2, , ] <- t(func_matrix)

  return(beta)
}


get_lambda <- function(col_name) {
  as.numeric(sub("aligned_lambda_", "", col_name))
}


fdasrvf_to_matrix <- function(aligned_curves, grid) {
  # Extract aligned curves from fdasrvf output (betan)
  beta_aligned <- aligned_curves$betan

  # Extract y-coordinates (dimension 2)
  f <- beta_aligned[2, , ]

  #TODO: fdasrvf seems to standardize to 0 mean
  # figure out how to add back the mean
  # centers <- aligned_curves$cent
  #
  # for(i in 1:ncol(f)) {
  #   f[, i] <- f[, i] + centers[2, i]
  # }

  # Return M x N matrix (transpose so time points are rows)
  t(f)
}


align_functions <- function(fun_curves, lambda = 0.0, parallel = FALSE,
                            max_iter = 10,
                            t_grid = seq(0, 1, length.out = 100)){

  # data format needed for fdasrvf library
  curves_arr <- tf_to_fdasrvf(fun_curves, curves, t_grid)

  # TODO: test which arguments lead to crashes and whether parallelizing
  # here is efficient
  data_aligned <- curve_srvf_align(curves_arr, mode = "O", rotated = TRUE,
                                   maxit = max_iter, parallel = parallel,
                                   lambda = lambda,
                                   scale = FALSE)

  # Then convert back to tf format
  result <- fdasrvf_to_matrix(data_aligned, t_grid)
  aligned_curves <- tfd(result, arg = t_grid)

  return(aligned_curves)
}


# Computes mean function, locations of peaks and valleys, and height and
# curvature of peaks
analyze_mean_function <- function(curves) {

  curves_mean <- mean(curves)
  #TODO: vary smoothing by arguments to function
  #curves_mean <- tf_smooth(curves_mean)
  slope <- tf_derive(curves_mean, order = 1)
  second_derivative <- tf_derive(curves_mean, order = 2)

  # Find peaks (local maxima) - where slope changes from positive to negative
  peaks <- tf_where(
    slope,
    value < 0 & dplyr::lag(value, 1, value[1]) > 0)

  # Find valleys (local minima) - where slope changes from negative to positive
  valleys <- tf_where(
    slope,
    value > 0 & dplyr::lag(value, 1, value[1]) < 0)

  # Compute height and negative curvature at peak locations
  # Accounting for reduction of arg at edges of interval when calculating
  # derivative
  peak_loc <- peaks[[1]]
  peak_loc <- peak_loc[peak_loc < max(tf_arg(second_derivative)) &
                         peak_loc > min(tf_arg(second_derivative))]
  peak_heights <- curves_mean[, peak_loc]
  peak_curvatures <- - second_derivative[, peak_loc]

  # ensuring non-negative curvatures and normalizing
  # is there a way to do this in tidyfun?
  peak_curvatures <- ifelse(peak_curvatures < 0, 0, peak_curvatures)
  max_curvature <- max(peak_curvatures)

  if (max_curvature > 0){
    peak_curvatures <- peak_curvatures / max_curvature
  }

  return(list(
    peaks = peaks[[1]],
    valleys = valleys[[1]],
    function_mean = curves_mean,
    peak_heights = peak_heights,
    peak_curvatures = peak_curvatures
  ))
}


# Finds which labeled peaks are in same location as peak in previous curve
peak_successor <- function(peaks1, valleys1, peaks2, valleys2, labels1) {
  # Helper function to compute ranges for peaks
  compute_peak_ranges <- function(peaks, valleys) {
    # Ensure valleys include endpoints
    all_valleys <- sort(unique(c(t_grid[0], valleys, t_grid[-1])))

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
  assign_labels_to_peaks <- function(peaks2, ranges, peaks1, labels1) {
    labels <- numeric(length(peaks2))

    for(i in seq_along(peaks2)) {
      # Find which ranges contain this peak
      in_range <- peaks2[i] >= ranges[,1] & peaks2[i] <= ranges[,2]
      matching_ranges <- which(in_range)

      if(length(matching_ranges) > 0) {
        if(length(matching_ranges) > 1) {
          # If multiple ranges match, choose the closest peak
          distances <- abs(peaks1[matching_ranges] - peaks2[i])
          closest_idx <- matching_ranges[which.min(distances)]
          labels[i] <- labels1[closest_idx]
        } else {
          labels[i] <- labels1[matching_ranges]
        }
      }
    }
    labels
  }

  # Helper function to resolve overlapping labels
  resolve_overlapping_labels <- function(labels, peaks1, peaks2, labels1) {
    unique_labels <- unique(labels[labels > 0])

    for(label in unique_labels) {
      duplicates <- which(labels == label)
      if(length(duplicates) > 1) {
        # Keep the closest peak and reset others
        orig_peak <- peaks1[which(labels1 == label)]
        distances <- abs(orig_peak - peaks2[duplicates])
        keep_idx <- duplicates[which.min(distances)]
        remove_idx <- duplicates[duplicates != keep_idx]
        labels[remove_idx] <- 0
      }
    }
    labels
  }

  # Main function logic
  if(length(peaks1) == 0) {
    # If no peaks in first function, assign new labels to all peaks in second
    labelMax <- max(labels1, 0)
    labels2 <- labelMax + seq_along(peaks2)
    return(list(
      labels = labels2,
      labelMax = labelMax + length(peaks2)
    ))
  }

  # Compute ranges for peaks in first function
  ranges <- compute_peak_ranges(peaks1, valleys1)

  # Assign initial labels to peaks in second function
  labels2 <- assign_labels_to_peaks(peaks2, ranges, peaks1, labels1)

  # Resolve any overlapping labels
  labels2 <- resolve_overlapping_labels(labels2, peaks1, peaks2, labels1)

  # Assign new labels to unmatched peaks
  labelMax <- max(labels1)
  unmatched <- which(labels2 == 0)
  if(length(unmatched) > 0) {
    labels2[unmatched] <- labelMax + seq_along(unmatched)
    labelMax <- labelMax + length(unmatched)
  }

  # Return both new labels and updated labelMax
  return(list(
    labels = labels2,
    labelMax = labelMax
  ))
}


# Finds which peaks are significant (curvature above threshold)
get_significant_peaks <- function(curvatures, labels, tau){
  sig_curvatures <- unlist(curvatures > tau)
  sig_labels <- labels[sig_curvatures]
  sig_labels
}


# Finds which peaks are significant for a large enough interval
find_persistent_peaks <- function(significant_peaks, ratio_threshold) {
  # Create a map of all unique labels and their appearances
  all_labels <- unique(unlist(significant_peaks))

  # Get lambda values from names
  lambda_values <- unlist(lapply(names(significant_peaks), function(name){
    get_lambda(name)
  }))

  # Initialize list to store intervals for each label
  label_intervals <- list()

  # For each label, find all lambda values where it appears
  for(label in all_labels) {
    # Find in which lambda indices this label appears
    appearances <- which(sapply(significant_peaks, function(peaks) label %in% peaks))

    if(length(appearances) > 0) {
      label_intervals[[as.character(label)]] <- lambda_values[appearances]
    }
  }

  # Find the label with the most appearances
  max_appearances <- max(sapply(label_intervals, length))

  # Find persistent peaks (those appearing in >= ratio_threshold * max_appearances lambdas)
  persistent_labels <- names(which(sapply(label_intervals, length) >= ratio_threshold * max_appearances))
  persistent_labels <- as.numeric(persistent_labels)

  return(list(
    persistent_labels = persistent_labels,
    intervals = label_intervals
  ))
}


# main function
peak_persistance_diagram <- function(curves, t_grid, max_lambda = 2, n_lambda = 10,
                                     sig_threshold = 0.03, pers_threshold = 0.28){
  # Create dataframe with function curves
  curves <- curve_data$curves
  curves <- tfb(curves, basis = "spline", bs = "bs")
  curves_df <- data.frame(curves)
  colnames(curves_df) <- c("curves")

  # Grid of lambda values
  lambda_values <- seq(0, max_lambda, length.out = n_lambda)

  # Create column names that incorporate lambda values
  # Format: "aligned_lambda_X" where X is the lambda value
  col_names <- paste0("aligned_lambda_", formatC(lambda_values, format="f", digits=3))

  # Apply align_functions for each lambda and combine results
  aligned_list <- lapply(lambda_values, function(lambda) {
    print(paste("Alignment for:", lambda))
    aligned <- align_functions(curves_df, lambda = lambda, parallel = TRUE,
                               max_iter = 5, t_grid = curve_data$t_grid)
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
      labels1 = labels1
    )

    # Store new labels
    all_labels[[curr_name]] <- result$labels
  }

  # Calculate which peaks are significant
  significant_peaks <- lapply(col_names, function(name){
    get_significant_peaks(peak_curvatures[[name]], all_labels[[name]], sig_threshold)
  })

  names(significant_peaks) <- col_names


  # Calculate persistance of peaks
  res <- find_persistent_peaks(significant_peaks, pers_threshold)

  # Persistance diagram in barchart form
  persistence_plot_bc <- create_persistence_diagram_bc(res$intervals)

  return(list(
    bc = persistence_plot_bc,
    significant_peaks = significant_peaks,
    intervals = res$intervals
  ))
}


# Peak-persistence diagram in barchart form
create_persistence_diagram_bc <- function(intervals) {
  # Create a dataframe with start and end points for each peak's persistence
  plot_data <- lapply(names(intervals), function(label) {
    lambdas <- sort(intervals[[label]])
    if(length(lambdas) > 0) {
      data.frame(
        peak_label = as.numeric(label),
        lambda_start = min(lambdas),
        lambda_end = max(lambdas)
      )
    }
  }) %>% bind_rows()

  # Calculate the width needed for bars to touch
  # We'll use half the minimum difference between peak labels
  peak_labels <- sort(unique(plot_data$peak_label))
  if(length(peak_labels) > 1) {
    bar_width <- min(diff(peak_labels))/2
  } else {
    bar_width <- 0.5  # default width if only one peak
  }

  # Create the plot with vertical lines showing persistence intervals
  ggplot(plot_data) +
    # Use linewidth parameter to create touching bars
    geom_segment(aes(x = peak_label, xend = peak_label,
                     y = lambda_start, yend = lambda_end),
                 linewidth = bar_width * 20, color = "blue") +
    # Set x-axis to show only the actual peak labels
    scale_x_continuous(breaks = peak_labels,
                       labels = as.integer(peak_labels)) +
    theme_minimal() +
    labs(
      title = "Peak Persistence Diagram",
      subtitle = "Vertical bars show lambda ranges where peaks remain significant",
      x = "Peak Label",
      y = "Lambda Value"
    ) +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
}
