library(tidyfun)
library(ggplot2)
library(tidyverse)


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


# Helper function to create location matrices
create_location_matrix <- function(locs, labels, lambda_values) {
  label_max = max(unlist(labels))
  n_lams <- length(lambda_values)

  # Initialize matrices
  location_matrix <- matrix(NA, nrow = n_lams, ncol = label_max)

  # Fill location matrix
  for (i in 1:n_lams) {
    if (length(labels[[i]]) > 0) {
      location_matrix[i, labels[[i]]] <- locs[[i]]
    }
  }

  return(location_matrix)
}


# Main function to draw PPD surface
create_persistance_diagram_sf <- function(t_grid, lambda_values, mean_functions,
                                          peak_locs, labels, sig_peaks) {

  location_matrix <- create_location_matrix(peak_locs, labels, lambda_values)

  # Each mean function is a tf with one function
  mfn_list <- lapply(mean_functions, function(item){
    item[[1]]
  })

  # warning thrown here, not sure how to remove
  mfn_df <- tfd(mfn_list, arg = t_grid)

  plot_data <- as.data.frame(mfn_df, unnest = TRUE, arg = t_grid)
  plot_data$lambda <- get_lambda(plot_data$id)

  # Create the plot
  plot <- ggplot(plot_data, aes(x = arg, y = lambda)) +
    geom_tile(aes(fill = value)) +
    # Use viridis color scale (similar to MATLAB's parula)
    scale_fill_viridis_c(
      name = "g̃λ(t)",
      option = "viridis"
    ) +
    # Adjust the aspect ratio to make it look more like MATLAB's output
    coord_fixed(ratio = diff(range(t_grid)) / diff(range(lambda_values))) +
    # Add labels and theme
    labs(
      x = "t",
      y = "λ"
    ) +
    theme_minimal() +
    # Fine-tune the appearance
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )

  # Create a data frame for peak locations
  # We'll process the indicator matrix to get x-y coordinates for each peak
  peak_data <- data.frame()
  label_data <- data.frame()

  # For each column (peak) in the indicator matrix
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
        group = paste0("peak_", peak_idx),
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
      geom_path(data = peak_data,
                aes(x = x, y = y, group = group),
                color = "grey50",
                linewidth = 2) +
      # Points marking peak positions
      geom_point(data = peak_data,
                 aes(x = x, y = y),
                 color = ifelse(peak_data$significant, "black", "grey50"),
                 size = 2) +
      # Add labels at the end of each peak line
      geom_text(data = label_data,
                aes(x = x, y = y, label = label),
                color = "grey20",
                size = 4,
                nudge_x = diff(range(t_grid))/50,  # Slight offset to the right
                nudge_y = diff(range(lambda_values))/50,  # Slight offset upward
                fontface = "bold")
  }

  return(plot)
}
