library(tidyfun)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(grid)


# Plotting a set of function curves, optionally with base function g in red
plot_simulated_curves <- function(curves, t_grid, g = NULL) {

  curve_data <- as.data.frame(curves)
  colnames(curve_data) <- c("curves")

  # plotting with function in red if base function has been passed as argument
  if (!is.null(g)){
    base_curve <- g(t_grid)
    curve_data$base_curve <- tfd(base_curve, arg = t_grid)

    curves_plot <- ggplot(curve_data, aes(y = curves)) +
      geom_spaghetti(alpha = 0.3) +
      geom_spaghetti(aes(y = base_curve), color = "red", linewidth = 1.25) +
      theme_minimal() +
      labs(title = "Functional Curves with Noise",
           x = "Time", y = "Value")

  } else {
    curves_plot <- ggplot(curve_data, aes(y = curves)) +
      geom_spaghetti(alpha = 0.3) +
      theme_minimal() +
      labs(title = "Functional Curves with Noise",
           x = "Time", y = "Value")
  }

  # TODO maybe: different colors for curves?

  # Alternative plot with base R
  # plot(curve_data, type = "spaghetti")
  # lines(base_curve, col = "red")

  curves_plot

}


# Plotting the base function with only warping function, amplitude factors or
# noise function applied, visualizing the impact of each
visualize_parameter_impacts_tf <- function(functional_data) {
  t_grid <- functional_data$t_grid
  g <- functional_data$base_function
  base_curve <- g(t_grid)

  curve_data <- as.data.frame(tfd(functional_data$curves, t_grid))
  colnames(curve_data) <- c("curves")

  # Calculate the curves with only one of the noise parameters/function each applies
  curve_data$base_curve <- tfd(base_curve, t_grid)
  curve_data$a_curves <- curve_data$base_curve * functional_data$amplitude_params
  curve_data$warped_curves <- tfd(g(functional_data$grid_warped), t_grid)
  curve_data$noise_curves <- tfd(functional_data$noise_grid) + curve_data$base_curve

  # Create individual plots
  plot_a <- ggplot(curve_data, aes(y = a_curves)) +
    geom_spaghetti(alpha = 0.5) +
    geom_spaghetti(aes(y = base_curve), color = "red", linewidth = 1.25) +
    theme_minimal() +
    labs(title = "Multiplicative Factors (a_i)",
         x = "Time", y = "Value")

  # Create individual plots
  plot_warped <- ggplot(curve_data, aes(y = warped_curves)) +
    geom_spaghetti(alpha = 0.5) +
    geom_spaghetti(aes(y = base_curve), color = "red", linewidth = 1.25) +
    theme_minimal() +
    labs(title = "Time Warping (γ_i)",
         x = "Time", y = "Value")

  plot_noise <- ggplot(curve_data, aes(y = noise_curves)) +
    geom_spaghetti(alpha = 0.5) +
    geom_spaghetti(aes(y = base_curve), color = "red", linewidth = 1.25) +
    theme_minimal() +
    labs(title = "Noise (ε_i)",
         x = "Time", y = "Value")

  # Combine plots
  grid.arrange(plot_a, plot_warped, plot_noise, ncol = 3,
               top = textGrob("Impact of parameters",
                              gp = gpar(fontsize = 16, fontface = "bold")))

}


# Plotting a set of time warping functions
visualize_warping_functions <- function(warp_functions){

  warp_functions_df <- data.frame(warp_functions)
  colnames(warp_functions_df) <- c("funcs")

  ggplot(warp_functions_df, aes(y = funcs)) +
    geom_spaghetti() +
    geom_abline(linetype = "dotted",
                col = "red",
                linewidth = 1.25)

  # TODO: line through origin should only be defined between 0 and 1

}
