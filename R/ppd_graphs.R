#' Create Peak Persistence Diagram in Barchart Form
#'
#' @param intervals List where names are peak labels and values are vectors of lambda values
#'                  where the peak is significant
#'
#' @return A ggplot2 object representing the persistence diagram
#'
#' @import ggplot2
#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Assuming 'peak_intervals' is a list of significant lambda values for each peak
#' persistence_plot <- create_persistence_diagram_bc(peak_intervals)
#' print(persistence_plot)
#' }
create_persistence_diagram_bc <- function(intervals) {
  # Input validation
  if (missing(intervals)) {
    stop("'intervals' is required")
  }
  if (!is.list(intervals)) {
    stop("'intervals' must be a list")
  }
  if (length(intervals) == 0) {
    stop("'intervals' cannot be empty")
  }
  if (is.null(names(intervals)) || any(names(intervals) == "")) {
    stop("'intervals' must be a named list with peak labels as names")
  }

  # Check if dplyr is available
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is needed for this function to work. Please install it.")
  }

  # Check if ggplot2 is available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.")
  }

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
  }) %>% dplyr::bind_rows()

  # Check if we have valid plot data
  if (nrow(plot_data) == 0) {
    stop("No valid intervals found in the input data")
  }

  # Calculate the width needed for bars to touch
  # We'll use half the minimum difference between peak labels
  peak_labels <- sort(unique(plot_data$peak_label))
  if(length(peak_labels) > 1) {
    bar_width <- min(diff(peak_labels))/2
  } else {
    bar_width <- 0.5  # default width if only one peak
  }

  # Create the plot with vertical lines showing persistence intervals
  ggplot2::ggplot(plot_data) +
    # Use linewidth parameter to create touching bars
    ggplot2::geom_segment(ggplot2::aes(x = .data$peak_label, xend = .data$peak_label,
                                       y = .data$lambda_start, yend = .data$lambda_end),
                          linewidth = bar_width * 20, color = "blue") +
    # Set x-axis to show only the actual peak labels
    ggplot2::scale_x_continuous(breaks = peak_labels,
                                labels = as.integer(peak_labels)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Peak Persistence Diagram",
      subtitle = "Vertical bars show lambda ranges where peaks remain significant",
      x = "Peak Label",
      y = "Lambda Value"
    ) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_line(color = "grey95"),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10)
    )
}

#' Create Location Matrix for Peak Persistence Diagram
#'
#' @param locs List of peak locations for each lambda value
#' @param labels List of peak labels for each lambda value
#' @param lambda_values Numeric vector of lambda values
#'
#' @return A matrix where rows correspond to lambda values and columns to peak labels
#'
#' @keywords internal
create_location_matrix <- function(locs, labels, lambda_values) {
  # Input validation
  if (missing(locs)) {
    stop("'locs' is required")
  }
  if (missing(labels)) {
    stop("'labels' is required")
  }
  if (missing(lambda_values)) {
    stop("'lambda_values' is required")
  }

  # Type checking
  if (!is.list(locs)) {
    stop("'locs' must be a list")
  }
  if (!is.list(labels)) {
    stop("'labels' must be a list")
  }
  if (!is.numeric(lambda_values)) {
    stop("'lambda_values' must be a numeric vector")
  }

  # Length checking
  n_lams <- length(lambda_values)
  if (length(locs) != n_lams) {
    stop("'locs' must have the same length as 'lambda_values'")
  }
  if (length(labels) != n_lams) {
    stop("'labels' must have the same length as 'lambda_values'")
  }

  # Check if we have any labels
  all_labels_flat <- unlist(labels)
  if (length(all_labels_flat) == 0) {
    return(matrix(NA, nrow = n_lams, ncol = 0))
  }

  # Find maximum label
  label_max <- max(all_labels_flat, na.rm = TRUE)

  # Initialize matrices
  location_matrix <- matrix(NA, nrow = n_lams, ncol = label_max)

  # Fill location matrix
  for (i in 1:n_lams) {
    if (length(labels[[i]]) > 0) {
      # Validate that labels and locations match
      if (length(labels[[i]]) != length(locs[[i]])) {
        stop(paste("Mismatch between labels and locations length for lambda index", i))
      }

      # Check for values outside valid range
      if (any(labels[[i]] > label_max)) {
        stop(paste("Label value exceeds maximum label for lambda index", i))
      }
      if (any(labels[[i]] <= 0)) {
        stop(paste("Label value must be positive for lambda index", i))
      }

      # Assign locations to matrix
      location_matrix[i, labels[[i]]] <- locs[[i]]
    }
  }

  return(location_matrix)
}


#' Create Peak Persistence Diagram Surface Plot
#'
#' @param t_grid Numeric vector representing the time grid
#' @param lambda_values Numeric vector of lambda values
#' @param mean_functions List of mean functions for each lambda value
#' @param peak_locs List of peak locations for each lambda value
#' @param labels List of peak labels for each lambda value
#' @param sig_peaks List of significant peak indices for each lambda value
#'
#' @return A ggplot2 object representing the persistence diagram surface
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @import tf
#'
#' @examples
#' \dontrun{
#' # Assuming all input parameters are available
#' surface_plot <- create_persistance_diagram_sf(t_grid, lambda_values,
#'                                              mean_functions, peak_locs,
#'                                              labels, sig_peaks)
#' print(surface_plot)
#' }
create_persistance_diagram_sf <- function(t_grid, lambda_values, mean_functions,
                                          peak_locs, labels, sig_peaks) {
  # Input validation
  if (missing(t_grid)) {
    stop("'t_grid' is required")
  }
  if (missing(lambda_values)) {
    stop("'lambda_values' is required")
  }
  if (missing(mean_functions)) {
    stop("'mean_functions' is required")
  }
  if (missing(peak_locs)) {
    stop("'peak_locs' is required")
  }
  if (missing(labels)) {
    stop("'labels' is required")
  }
  if (missing(sig_peaks)) {
    stop("'sig_peaks' is required")
  }
  # Type checking
  if (!is.numeric(t_grid)) {
    stop("'t_grid' must be a numeric vector")
  }
  if (!is.numeric(lambda_values)) {
    stop("'lambda_values' must be a numeric vector")
  }
  if (!is.list(mean_functions)) {
    stop("'mean_functions' must be a list")
  }
  if (!is.list(peak_locs)) {
    stop("'peak_locs' must be a list")
  }
  if (!is.list(labels)) {
    stop("'labels' must be a list")
  }
  if (!is.list(sig_peaks)) {
    stop("'sig_peaks' must be a list")
  }
  # Length checking
  n_lams <- length(lambda_values)
  if (length(mean_functions) != n_lams) {
    stop("'mean_functions' must have the same length as 'lambda_values'")
  }
  if (length(peak_locs) != n_lams) {
    stop("'peak_locs' must have the same length as 'lambda_values'")
  }
  if (length(labels) != n_lams) {
    stop("'labels' must have the same length as 'lambda_values'")
  }
  if (length(sig_peaks) != n_lams) {
    stop("'sig_peaks' must have the same length as 'lambda_values'")
  }
  # Check function existence
  if (!exists("get_lambda")) {
    stop("Function 'get_lambda' is required but not available")
  }
  # Check packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.")
  }

  # Create location matrix
  tryCatch({
    location_matrix <- create_location_matrix(peak_locs, labels, lambda_values)
  }, error = function(e) {
    stop("Error creating location matrix: ", e$message)
  })
  # Each mean function is a tf with one function
  mfn_list <- tryCatch({
    lapply(mean_functions, function(item) {
      item[[1]]
    })
  }, error = function(e) {
    stop("Error processing mean functions: ", e$message)
  })
  # warning thrown here, not sure how to remove
  tryCatch({
    mfn_df <- tfd(mfn_list, arg = t_grid)
  }, error = function(e) {
    stop("Error converting mean functions to tfd format: ", e$message)
  })

  plot_data <- tryCatch({
    as.data.frame(mfn_df, unnest = TRUE, arg = t_grid)
  }, error = function(e) {
    stop("Error creating plot data from mean functions: ", e$message)
  })

  plot_data$lambda <- get_lambda(plot_data$id)
  # Create the plot
  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$arg,
                                                  y = .data$lambda)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$value)) +
    # Use viridis color scale (similar to MATLAB's parula)
    ggplot2::scale_fill_viridis_c(
      name = "g(t)",
      option = "viridis"
    ) +
    # Adjust the aspect ratio to make it look more like MATLAB's output
    ggplot2::coord_fixed(ratio = diff(range(t_grid)) / diff(range(lambda_values))) +
    # Add labels and theme
    ggplot2::labs(
      x = "t",
      y = "Lambda"
    ) +
    ggplot2::theme_minimal() +
    # Fine-tune the appearance
    ggplot2::theme(
      legend.position = "right",
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 10)
    )
  # Create a data frame for peak locations
  # We'll process the indicator matrix to get x-y coordinates for each peak
  peak_data <- data.frame()
  label_data <- data.frame()
  # For each column (peak) in the location matrix
  for (peak_idx in 1:ncol(location_matrix)) {
    # Find which lambda values (rows) have this peak
    lambda_indices <- which(!is.na(location_matrix[, peak_idx]))
    if (length(lambda_indices) > 0) {
      # Get x and y coordinates for this peak's path
      x_coords <- location_matrix[lambda_indices, peak_idx]
      y_coords <- lambda_values[lambda_indices]
      # Determine significance
      # For each point along this peak's path, determine if it's significant
      # by checking the significant_peaks list for each lambda value
      is_significant <- sapply(lambda_indices, function(lambda_idx) {
        # Get the significant peaks for this lambda value
        sig_peaks_at_lambda <- sig_peaks[[lambda_idx]]
        # Check if this peak is in the significant peaks for this lambda
        peak_idx %in% sig_peaks_at_lambda
      })
      # Create temporary data frame for this peak's path
      temp_df <- data.frame(
        # x coordinates come from t_grid positions indicated in the matrix
        x = x_coords,
        # y coordinates are the corresponding lambda values
        y = y_coords,
        # Group identifier for connecting points
        group = paste0("peak", peak_idx),
        # mask indicating whether peak is significant
        significant = is_significant
      )
      peak_data <- rbind(peak_data, temp_df)
      # Create label data frame for the end of this peak's path
      # We'll use the last point (highest lambda value) for the label
      # Label color based on significance of final lambda
      label_df <- data.frame(
        x = x_coords[length(x_coords)],
        y = y_coords[length(y_coords)],
        label = as.character(peak_idx),
        significant = is_significant[length(is_significant)]
      )
      label_data <- rbind(label_data, label_df)
    }
  }
  # Add peak markers and connections if we have peak data
  if (nrow(peak_data) > 0) {
    plot <- plot +
      # Lines connecting peak positions
      ggplot2::geom_path(data = peak_data,
                         ggplot2::aes(x = .data$x, y = .data$y, group = .data$group),
                         color = "grey50",
                         linewidth = 2) +
      # Points marking peak positions
      ggplot2::geom_point(data = peak_data,
                          ggplot2::aes(x = .data$x, y = .data$y),
                          color = ifelse(peak_data$significant, "black", "grey50"),
                          size = 2) +
      # Add labels at the end of each peak line
      ggplot2::geom_text(data = label_data,
                         ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
                         color = "grey20",
                         size = 4,
                         nudge_x = diff(range(t_grid))/50,  # Slight offset to the right
                         nudge_y = diff(range(lambda_values))/50,  # Slight offset upward
                         fontface = "bold")
  }
  return(plot)
}
