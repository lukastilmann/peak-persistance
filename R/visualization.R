#' Plot Simulated Functional Curves
#'
#' Creates a visualization of simulated functional curves. If a base function
#' is provided, it will be plotted in red to show the original function from
#' which the curves were derived.
#'
#' @param curves A tidyfun data object containing functional curves.
#' @param t_grid Numeric vector of grid points on which the curves are evaluated.
#' @param g Function or NULL. The base function used to generate the curves.
#'        If provided, it will be plotted in red. Default is NULL.
#'
#' @return A ggplot object displaying the curves.
#' @export
#'
#' @examples
#' # Create a simple base function
#' g <- function(x) sin(2*pi*x)
#'
#' # Generate sample curves
#' t_grid <- seq(0, 1, length.out = 100)
#' curves_data <- generate_functional_curves(
#'   n = 20,
#'   num_points = 100,
#'   g = g,
#'   seed = 123
#' )
#'
#' # Plot the curves with base function
#' plot_simulated_curves(curves_data$curves, curves_data$t_grid, g)
plot_simulated_curves <- function(curves, t_grid, g = NULL) {
  # Input validation
  if (!inherits(curves, "tf")) {
    stop("'curves' must be a tidyfun data object (tf)")
  }

  if (!is.numeric(t_grid)) {
    stop("'t_grid' must be a numeric vector")
  }

  if (!is.null(g) && !is.function(g)) {
    stop("'g' must be a function or NULL")
  }

  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }

  if (!requireNamespace("tidyfun", quietly = TRUE)) {
    stop("Package 'tidyfun' is required. Please install it.")
  }

  # Convert to data frame
  curve_data <- as.data.frame(curves)
  colnames(curve_data) <- c("curves")

  # Plotting with base function in red if provided
  if (!is.null(g)) {
    # Evaluate base function on grid
    base_curve <- g(t_grid)

    # Convert to tidyfun object for consistent plotting
    curve_data$base_curve <- tf::tfd(base_curve, arg = t_grid)

    # Create plot with base function
    curves_plot <- ggplot2::ggplot(curve_data, ggplot2::aes(y = curves)) +
      tidyfun::geom_spaghetti(alpha = 0.3) +
      tidyfun::geom_spaghetti(ggplot2::aes(y = base_curve), color = "red", linewidth = 1.25) +
      #ggplot2::theme_minimal() +
      ggplot2::labs(title = "Functional Curves with Noise",
                    x = "Time", y = "Value")
  } else {
    # Create plot without base function
    curves_plot <- ggplot2::ggplot(curve_data, ggplot2::aes(y = curves)) +
      tidyfun::geom_spaghetti(alpha = 0.3) +
      #ggplot2::theme_minimal() +
      ggplot2::labs(title = "Functional Curves with Noise",
                    x = "Time", y = "Value")
  }

  return(curves_plot)
}

#' Visualize Parameter Impacts on Functional Curves
#'
#' Creates a multi-panel plot that shows the individual effects of amplitude factors,
#' time warping, and noise on the base function. This visualization helps understand
#' how each component contributes to the variability in the simulated curves.
#'
#' @param functional_data A list containing the output from generate_functional_curves(),
#'        including curves, t_grid, base_function, amplitude_params, grid_warped, and
#'        noise_grid.
#'
#' @return A grid of ggplot objects displaying the parameter impacts.
#' @export
#'
#' @examples
#' # Create a simple base function
#' g <- function(x) sin(2*pi*x)
#'
#' # Generate sample curves with parameter variations
#' curves_data <- generate_functional_curves(
#'   n = 20,
#'   g = g,
#'   sigma_amplitude = 0.1,
#'   warping = "flexible",
#'   noise_to_signal = 0.05,
#'   seed = 123
#' )
#'
#' # Visualize parameter impacts
#' visualize_parameter_impacts_tf(curves_data)
visualize_parameter_impacts_tf <- function(functional_data) {
  # Input validation
  if (!is.list(functional_data)) {
    stop("'functional_data' must be a list")
  }

  required_elements <- c("t_grid", "base_function", "curves",
                         "amplitude_params", "grid_warped", "noise_grid")
  missing_elements <- setdiff(required_elements, names(functional_data))

  if (length(missing_elements) > 0) {
    stop("Missing required elements in 'functional_data': ",
         paste(missing_elements, collapse = ", "),
         ". Make sure it's the output from generate_functional_curves().")
  }

  if (!is.function(functional_data$base_function)) {
    stop("'functional_data$base_function' must be a function")
  }

  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }

  if (!requireNamespace("tidyfun", quietly = TRUE)) {
    stop("Package 'tidyfun' is required. Please install it.")
  }

  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Package 'gridExtra' is required. Please install it.")
  }

  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package 'grid' is required. Please install it.")
  }

  # Extract data from input
  t_grid <- functional_data$t_grid
  g <- functional_data$base_function

  # Evaluate base function on grid
  base_curve <- g(t_grid)

  # Convert to data frame with tidyfun
  curve_data <- as.data.frame(tf::tfd(functional_data$curves, t_grid))
  colnames(curve_data) <- c("curves")

  # Calculate curves with only one parameter/function applied at a time
  curve_data$base_curve <- tf::tfd(base_curve, t_grid)
  curve_data$a_curves <- curve_data$base_curve * functional_data$amplitude_params
  curve_data$warped_curves <- tf::tfd(g(functional_data$grid_warped), t_grid)
  curve_data$noise_curves <- tf::tfd(functional_data$noise_grid) + curve_data$base_curve

  # Create individual plots
  plot_a <- ggplot2::ggplot(curve_data, ggplot2::aes(y = .data$a_curves)) +
    tidyfun::geom_spaghetti(alpha = 0.5) +
    tidyfun::geom_spaghetti(ggplot2::aes(y = base_curve), color = "red", linewidth = 1.25) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Multiplicative Factors (a_i)",
                  x = "Time", y = "Value")

  plot_warped <- ggplot2::ggplot(curve_data, ggplot2::aes(y = .data$warped_curves)) +
    tidyfun::geom_spaghetti(alpha = 0.5) +
    tidyfun::geom_spaghetti(ggplot2::aes(y = base_curve), color = "red", linewidth = 1.25) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Time Warping",
                  x = "Time", y = "Value")

  plot_noise <- ggplot2::ggplot(curve_data, ggplot2::aes(y = .data$noise_curves)) +
    tidyfun::geom_spaghetti(alpha = 0.5) +
    tidyfun::geom_spaghetti(ggplot2::aes(y = base_curve), color = "red", linewidth = 1.25) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Noise",
                  x = "Time", y = "Value")

  # Combine plots
  gridExtra::grid.arrange(plot_a, plot_warped, plot_noise, ncol = 3,
                          top = grid::textGrob("Impact of parameters",
                                               gp = grid::gpar(fontsize = 16, fontface = "bold")))
}

#' Visualize Warping Functions
#'
#' Creates a plot of time warping functions to show how they distort the
#' time axis. The identity function (y = x) is shown as a dotted red line
#' for reference.
#'
#' @param warp_functions A tidyfun data object containing warping functions.
#'
#' @return A ggplot object displaying the warping functions.
#' @export
visualize_warping_functions <- function(warp_functions) {
  # Input validation
  if (!inherits(warp_functions, "tf")) {
    stop("'warp_functions' must be a tidyfun data object (tf)")
  }

  # Check for required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Please install it.")
  }

  if (!requireNamespace("tidyfun", quietly = TRUE)) {
    stop("Package 'tidyfun' is required. Please install it.")
  }

  # Convert to data frame
  warp_functions_df <- as.data.frame(warp_functions)
  colnames(warp_functions_df) <- c("funcs")

  # Create plot
  warp_plot <- ggplot2::ggplot(warp_functions_df, ggplot2::aes(y = .data$funcs)) +
    tidyfun::geom_spaghetti() +
    ggplot2::geom_abline(linetype = "dotted", col = "red", linewidth = 1.25) +
    ggplot2::xlim(0, 1) +     # Limit x-axis to [0,1]
    ggplot2::ylim(0, 1) +     # Limit y-axis to [0,1]
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Warping Functions",
                  x = "Original Time", y = "Warped Time")

  return(warp_plot)
}


