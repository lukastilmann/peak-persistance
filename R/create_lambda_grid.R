#' Create Parameter Grids with Different Spacings
#'
#' @description
#' Creates a grid of parameter values using different spacing methods,
#' typically for regularization parameters in statistical models.
#'
#' @param max_value Numeric. The maximum value in the grid.
#' @param n_points Integer. The number of points in the grid. Default is 10.
#' @param spacing Either "log" for logarithmic spacing or a numeric value for power-based spacing.
#'                If "log", points are spaced logarithmically between a small value and max_value.
#'                If numeric (e.g., 2), points are spaced according to value^(1/spacing). Default is "log".
#'
#' @return A numeric vector of length n_points containing the parameter grid with the specified spacing.
#'         The grid always includes 0 when spacing="log".
#'
#' @examples
#' # Create a logarithmically spaced grid
#' log_grid <- create_lambda_grid(10, n_points = 5, spacing = "log")
#'
#' # Create a square root spaced grid (spacing = 2)
#' sqrt_grid <- create_lambda_grid(10, n_points = 5, spacing = 2)
#'
#' @export
create_lambda_grid <- function(max_value = 2, n_points = 10,
                               lambda_grid_spacing = c("log", "sqrt")) {
  # Input validation
  checkmate::assert_number(max_value, lower = 0, finite = TRUE)
  lambda_grid_spacing <- match.arg(lambda_grid_spacing)

  if (lambda_grid_spacing == "log") {
    min_value <- 1e-4  # Small positive number close to zero
    return(c(0, exp(seq(log(min_value), log(max_value), length.out = n_points - 1))))
  } else {
    # For square root spacing
    return(seq(0, (max_value)^(1/2), length.out = n_points)^2)
  }
}


