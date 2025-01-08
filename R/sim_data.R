library(tidyfun)
library(ggplot2)
library(gtools)


# Flexible function generator for mathematical benchmark functions
# int indicates x^n as basis function, "normal" a normal density,
# and "sine" a sine curve
generate_benchmark_function <- function(bases,
                                        coefficients,
                                        shifts = NULL,
                                        scale_factors = NULL) {
  # Default to 0 for shifts and 1 for amplitudes if not provided
  if (is.null(shifts)) shifts <- rep(0, length(bases))
  if (is.null(scale_factors)) amplitudes <- rep(1, length(bases))

  # Validate input lengths
  if (length(bases) != length(coefficients) ||
      length(bases) != length(shifts) ||
      length(bases) != length(scale_factors)) {
    stop("bases, coefficients, shifts, and scale_factors must have the same length")
  }

  # Create the combined function
  function(x) {
    # Initialize total result
    total_result <- 0

    # Iterate through each basis
    for (i in seq_along(bases)) {
      # Polynomial basis (numeric power)
      if (is.numeric(bases[[i]])) {
        # Compute transformed x
        x_transformed <- (x - shifts[i]) * scale_factors[i]

        # Compute term with specified power and coefficient
        term <- coefficients[[i]] * x_transformed^bases[[i]]

        total_result <- total_result + term
      }
      # Normal density basis
      else if (bases[[i]] == "normal") {
        # Compute transformed x
        x_transformed <- (x - shifts[i]) * scale_factors[i]

        # Calculate scaled normal density
        term <- coefficients[[i]] * dnorm(x_transformed)

        total_result <- total_result + term
      }
      # Sine basis
      else if (bases[[i]] == "sine") {
        # Compute transformed x
        x_transformed <- (x - shifts[i]) * scale_factors[i]

        # Calculate scaled sine
        term <- coefficients[[i]] * sin(x_transformed)

        total_result <- total_result + term
      }
      else {
        stop("Unsupported basis type")
      }
    }

    return(total_result)
  }
}


# Function to generate boundary-preserving warping functions
generate_warping_function <- function(n_points = 5,
                                      lambda = 0,
                                      seed = NULL) {
  # Set seed for reproducibility if provided
  if (!is.null(seed)) set.seed(seed)

  # Sample n points between 0 and 1 and order them
  sample_points <- sort(runif(n_points))

  # Calculate cumulative increments
  # Add 0 and 1 as boundary points
  full_points <- c(0, sample_points, 1)
  interval_sizes <- diff(full_points)

  # Sample increments from Dirichlet distribution
  # Add 1 to each interval size to ensure positive parameters
  dirichlet_params <- interval_sizes * (1 + lambda)
  increments <- rdirichlet(1, dirichlet_params)[1,]

  # Calculate cumulative function values
  cumulative_values <- c(0, cumsum(increments))

  # Normalize to ensure endpoint is 1
  cumulative_values <- cumulative_values / cumulative_values[length(cumulative_values)]

  # Create warping function using approx()
  # TODO: function should probably be smoother
  function(x) {
    # Ensure x is within [0, 1]
    x <- pmax(0, pmin(1, x))

    # Use linear interpolation
    approx(full_points, cumulative_values, xout = x, method = "linear")$y
  }
}


# Generate a given number of functional curves based on a base function g
generate_functional_curves <- function(n = 50,
                                       num_points = 100,
                                       g,
                                       warping = "simple",
                                       warping_lambda = 0,
                                       warping_points = 5,
                                       sigma_amplitude = 0.01,
                                       scale_factor = 0.01,
                                       nugget = 0.01) {
  # number of points in grid
  t_grid <- seq(0, 1, length.out = num_points)

  # Generate amplitude factors from N(1, sigma_amplitude^2)
  a_params <- rnorm(n, mean = 1, sd = sigma_amplitude)

  if (warping == "simple"){
    # Generate time-warping parameters from U[-1,1]
    z_params <- runif(n, -1, 1)

    # Compute warped grids directly
    grid_warped <- sapply(z_params, function(zi) {
      t_grid + zi * t_grid * (1 - t_grid)
    })
  } else if (warping == "flexible"){
    # Draw random warping function with more flexibility
    grid_warped <- sapply(1:n, function(i){
      generate_warping_function(n_points = warping_points,
                                lambda = warping_lambda)(t_grid)
    })
  } else {
    stop("Warping type not specified correctly.")
  }

  grid_warped <- t(grid_warped)
  warped_values <- tfd(g(grid_warped), t_grid)

  noise <- tf_rgp(n = n,
                  arg = num_points,
                  cov = "squareexp",
                  scale = scale_factor,
                  nugget = nugget)

  curves <- warped_values * a_params + noise

  return(list(
    curves = curves,
    t_grid = t_grid,
    base_function = g,
    amplitude_params = a_params,
    grid_warped = grid_warped,
    noise_grid = noise
  ))

}


