% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_lambda_grid.R
\name{create_lambda_grid}
\alias{create_lambda_grid}
\title{Create Parameter Grids with Different Spacings}
\usage{
create_lambda_grid(
  max_value = 2,
  n_points = 10,
  lambda_grid_spacing = c("log", "sqrt", "cubicrt")
)
}
\arguments{
\item{max_value}{Numeric. The maximum value in the grid.}

\item{n_points}{Integer. The number of points in the grid. Default is 10.}

\item{lambda_grid_spacing}{Either "log" for logarithmic spacing or "sqrt" for squareroot
spacing or "cubicrt" for cubic root spacing.
If "log", points are spaced logarithmically between a small value and max_value.
If a root, points are spaced according to value^(1/2)or value^(1/3).
Default is "log".}
}
\value{
A numeric vector of length n_points containing the parameter grid with the specified spacing.
The grid always includes 0 when spacing="log".
}
\description{
Creates a grid of parameter values using different spacing methods,
typically for regularization parameters in statistical models.
}
\examples{
# Create a logarithmically spaced grid
log_grid <- create_lambda_grid(10, n_points = 5, lambda_grid_spacing = "log")

# Create a square root spaced grid
sqrt_grid <- create_lambda_grid(10, n_points = 5, lambda_grid_spacing = "sqrt")

}
