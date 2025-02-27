#' Generate a Flexible Mathematical Benchmark Function
#'
#' Creates a function by combining multiple basis functions with coefficients,
#' shifts, and scale factors. Supported basis types include polynomial functions
#' (specified by numeric power), normal density functions, and sine functions.
#'
#' @param bases List of basis function types. Use numeric values for polynomial
#'        powers (e.g., 1 for linear, 2 for quadratic), "normal" for normal
#'        density, or "sine" for sine function.
#' @param coefficients Numeric vector of coefficients for each basis function.
#' @param shifts Numeric vector of horizontal shifts for each basis function.
#'        Defaults to 0 for each basis if NULL.
#' @param scale_factors Numeric vector of scaling factors for each basis function.
#'        Defaults to 1 for each basis if NULL.
#'
#' @return A function that takes a numeric input and returns the combined value
#'         of all basis functions.
#' @export
#'
#' @examples
#' # Create a function with quadratic and sine components
#' f <- generate_benchmark_function(
#'   bases = list(2, "sine"),
#'   coefficients = c(1.5, 2),
#'   shifts = c(0, pi/4),
#'   scale_factors = c(1, 2)
#' )
#'
#' # Evaluate at points
#' x <- seq(-2, 2, by = 0.1)
#' y <- f(x)
#' plot(x, y, type = "l")
#' @param bases List of basis function types. Use numeric values for polynomial
#'        powers (e.g., 1 for linear, 2 for quadratic), "normal" for normal
#'        density, or "sine" for sine function.
#' @param coefficients Numeric vector or list of coefficients for each basis function.
#' @param shifts Numeric vector or list of horizontal shifts for each basis function.
#'        Defaults to 0 for each basis if NULL.
#' @param scale_factors Numeric vector or list of scaling factors for each basis function.
#'        Defaults to 1 for each basis if NULL.
#'
#' @return A function that takes a numeric input and returns the combined value
#'         of all basis functions.
#' @export
#'
#' @examples
#' # Create a function with quadratic and sine components
#' f <- generate_benchmark_function(
#'   bases = list(2, "sine"),
#'   coefficients = c(1.5, 2),
#'   shifts = c(0, pi/4),
#'   scale_factors = c(1, 2)
#' )
#'
#' # Evaluate at points
#' x <- seq(-2, 2, by = 0.1)
#' y <- f(x)
#' plot(x, y, type = "l")
generate_benchmark_function <- function(bases,
                                        coefficients,
                                        shifts = NULL,
                                        scale_factors = NULL) {
  # Convert list inputs to vectors if needed
  if (is.list(coefficients)) coefficients <- unlist(coefficients)

  # Input validation
  if (!is.list(bases)) {
    stop("'bases' must be a list")
  }

  if (!is.numeric(coefficients)) {
    stop("'coefficients' must be numeric (or a list that can be unlisted to numeric)")
  }

  # Default to 0 for shifts and 1 for scale_factors if not provided
  if (is.null(shifts)) shifts <- rep(0, length(bases))
  if (is.null(scale_factors)) scale_factors <- rep(1, length(bases))

  # Convert list inputs to vectors if needed (after defaults are set)
  if (is.list(shifts)) shifts <- unlist(shifts)
  if (is.list(scale_factors)) scale_factors <- unlist(scale_factors)

  # Validate input types
  if (!is.numeric(shifts)) stop("'shifts' must be numeric (or a list that can be unlisted to numeric)")
  if (!is.numeric(scale_factors)) stop("'scale_factors' must be numeric (or a list that can be unlisted to numeric)")

  # Validate input lengths
  if (length(bases) != length(coefficients)) {
    stop("'bases' and 'coefficients' must have the same length")
  }
  if (length(bases) != length(shifts)) {
    stop("'bases' and 'shifts' must have the same length")
  }
  if (length(bases) != length(scale_factors)) {
    stop("'bases' and 'scale_factors' must have the same length")
  }

  # Validate basis types
  valid_basis_types <- c("normal", "sine")
  for (basis in bases) {
    if (!is.numeric(basis) && !(basis %in% valid_basis_types)) {
      stop("Each basis must be either numeric or one of: ",
           paste(valid_basis_types, collapse = ", "))
    }
  }

  # Create the combined function
  function(x) {
    # Check input type
    if (!is.numeric(x)) {
      stop("Input to the generated function must be numeric")
    }

    # Initialize total result
    total_result <- 0

    # Iterate through each basis
    for (i in seq_along(bases)) {
      # Compute transformed x
      x_transformed <- (x - shifts[i]) * scale_factors[i]

      # Polynomial basis (numeric power)
      if (is.numeric(bases[[i]])) {
        # Compute term with specified power and coefficient
        term <- coefficients[i] * x_transformed^bases[[i]]
      }
      # Normal density basis
      else if (bases[[i]] == "normal") {
        # Calculate scaled normal density
        term <- coefficients[i] * stats::dnorm(x_transformed)
      }
      # Sine basis
      else if (bases[[i]] == "sine") {
        # Calculate scaled sine
        term <- coefficients[i] * sin(x_transformed)
      }
      else {
        stop("Unsupported basis type: ", bases[[i]])
      }

      total_result <- total_result + term
    }

    return(total_result)
  }
}


#' Generate a Boundary-Preserving Warping Function
#'
#' Creates a monotone increasing function that maps `[0,1]` to `[0,1]`, preserving the
#' endpoints. The function is built by sampling points within the unit interval and
#' creating a monotone cubic spline interpolation between them.
#'
#' @param n_points Integer. Number of internal points to sample within (0,1).
#'        More points allow for more complex warping. Default is 5.
#' @param gamma Numeric. Controls the variability of the warping function.
#'        Higher values increase variability. Default is 0.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#'
#' @return A function that takes numeric input(s) in `[0,1]` and returns warped values
#'         in `[0,1]`.
#' @export
#'
#' @examples
#' # Create a warping function with default parameters
#' warp_fn <- generate_warping_function(seed = 123)
#'
#' # Plot the warping function
#' x <- seq(0, 1, by = 0.01)
#' y <- warp_fn(x)
#' plot(x, y, type = "l", main = "Warping Function",
#'      xlab = "Original", ylab = "Warped")
#' abline(0, 1, lty = 2)  # Add identity line for reference
generate_warping_function <- function(n_points = 5,
                                      gamma = 0,
                                      seed = NULL) {
  # Input validation
  if (!is.numeric(n_points) || length(n_points) != 1 ||
      n_points < 1 || n_points != round(n_points)) {
    stop("'n_points' must be a positive integer")
  }

  if (!is.numeric(gamma) || length(gamma) != 1 || gamma < 0) {
    stop("'gamma' must be a non-negative numeric value")
  }

  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1 ||
                         seed != round(seed))) {
    stop("'seed' must be NULL or an integer")
  }

  # Set seed for reproducibility if provided
  if (!is.null(seed)) set.seed(seed)

  # Sample n points between 0 and 1 and order them
  sample_points <- sort(stats::runif(n_points))

  # Add 0 and 1 as boundary points
  full_points <- c(0, sample_points, 1)
  interval_sizes <- diff(full_points)

  # Sample increments from Dirichlet distribution
  # Add 1 to each interval size to ensure positive parameters
  dirichlet_params <- interval_sizes * (1 + gamma)

  # Check for gtools package
  if (!requireNamespace("gtools", quietly = TRUE)) {
    stop("Package 'gtools' is required for the Dirichlet distribution. Please install it.")
  }

  increments <- gtools::rdirichlet(1, dirichlet_params)[1,]

  # Calculate cumulative function values
  cumulative_values <- c(0, cumsum(increments))

  # Normalize to ensure endpoint is 1
  cumulative_values <- cumulative_values / cumulative_values[length(cumulative_values)]

  # Create warping function using splinefun with Hyman's monotone cubic interpolation
  function(x) {
    # Input validation for the returned function
    if (!is.numeric(x)) {
      stop("Input to warping function must be numeric")
    }

    # Ensure x is within [0, 1]
    x <- pmax(0, pmin(1, x))

    # Monotone cubic spline interpolation
    stats::splinefun(full_points, cumulative_values, method = "hyman")(x)
  }
}


#' Generate Functional Curves Based on a Benchmark Function
#'
#' Creates a set of functional curves by applying amplitude variations, time
#' warping, and random noise to a base function. This can be used to generate
#' synthetic functional data for testing and simulation purposes.
#'
#' @param n Integer. Number of curves to generate. Default is 50.
#' @param num_points Integer. Number of points in the grid for each curve. Default is 100.
#' @param g Function. Base function to be warped and perturbed.
#' @param warping Character. Type of warping: "simple" for polynomial warping or
#'        "flexible" for more complex warping. Default is "simple".
#' @param warping_gamma Numeric. Controls the variability when using flexible warping.
#'        Higher values increase variability. Default is 0.
#' @param warping_points Integer. Number of internal points for flexible warping.
#'        Default is 5.
#' @param sigma_amplitude Numeric. Standard deviation of amplitude variation.
#'        Default is 0.01.
#' @param scale_factor Numeric. Scale factor for warping (not currently used).
#'        Default is 0.01.
#' @param nugget Numeric. Base noise level (not currently used). Default is 0.01.
#' @param noise_to_signal Numeric. Ratio of noise variance to signal variance.
#'        Default is 0.01.
#' @param seed Integer or NULL. Random seed for reproducibility. Default is NULL.
#'
#' @return A list containing:
#'   \item{curves}{Matrix of generated functional curves.}
#'   \item{t_grid}{Vector of grid points (0 to 1).}
#'   \item{base_function}{The base function used.}
#'   \item{amplitude_params}{Vector of amplitude parameters for each curve.}
#'   \item{grid_warped}{Matrix of warped grid points.}
#'   \item{noise_grid}{Matrix of noise added to each curve.}
#' @export
#'
#' @examples
#' # Create a simple benchmark function
#' base_fn <- generate_benchmark_function(
#'   bases = list(2, "sine"),
#'   coefficients = c(1, 0.5),
#'   shifts = c(0.5, 0),
#'   scale_factors = c(1, 6.28)
#' )
#'
#' # Generate 20 curves with simple warping
#' curves_data <- generate_functional_curves(
#'   n = 20,
#'   num_points = 100,
#'   g = base_fn,
#'   warping = "simple",
#'   sigma_amplitude = 0.05,
#'   noise_to_signal = 0.02,
#'   seed = 123
#' )
#'
#' # Plot the first 5 curves
#' plot(curves_data$curves, curves_data$t_grid)
generate_functional_curves <- function(g,
                                       n = 50,
                                       num_points = 100,
                                       warping = "simple",
                                       warping_gamma = 0,
                                       warping_points = 5,
                                       sigma_amplitude = 0.01,
                                       scale_factor = 0.01,
                                       nugget = 0.01,
                                       noise_to_signal = 0.01,
                                       seed = NULL) {
  # Input validation
  if (!is.function(g)) {
    stop("'g' must be a function")
  }

  if (!is.numeric(n) || length(n) != 1 || n < 1 || n != round(n)) {
    stop("'n' must be a positive integer")
  }

  if (!is.numeric(num_points) || length(num_points) != 1 ||
      num_points < 2 || num_points != round(num_points)) {
    stop("'num_points' must be an integer greater than 1")
  }

  if (!is.character(warping) || length(warping) != 1 ||
      !(warping %in% c("simple", "flexible"))) {
    stop("'warping' must be either 'simple' or 'flexible'")
  }

  if (!is.numeric(warping_gamma) || length(warping_gamma) != 1 || warping_gamma < 0) {
    stop("'warping_gamma' must be a non-negative numeric value")
  }

  if (!is.numeric(warping_points) || length(warping_points) != 1 ||
      warping_points < 1 || warping_points != round(warping_points)) {
    stop("'warping_points' must be a positive integer")
  }

  if (!is.numeric(sigma_amplitude) || length(sigma_amplitude) != 1 || sigma_amplitude < 0) {
    stop("'sigma_amplitude' must be a non-negative numeric value")
  }

  if (!is.numeric(noise_to_signal) || length(noise_to_signal) != 1 || noise_to_signal < 0) {
    stop("'noise_to_signal' must be a non-negative numeric value")
  }

  # Check for required packages
  if (!requireNamespace("tf", quietly = TRUE)) {
    stop("Package 'tf' is required. Please install it.")
  }

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required. Please install it.")
  }

  # Set seed for reproducibility if provided
  if (!is.null(seed)) set.seed(seed)

  # Create grid of points from 0 to 1
  t_grid <- seq(0, 1, length.out = num_points)

  # Generate amplitude factors from N(1, sigma_amplitude^2)
  a_params <- stats::rnorm(n, mean = 1, sd = sigma_amplitude)

  # Apply time warping based on specified method
  if (warping == "simple") {
    # Generate time-warping parameters from U[-1,1]
    z_params <- stats::runif(n, -1, 1)

    # Compute warped grids using polynomial warping
    grid_warped <- sapply(z_params, function(zi) {
      t_grid + zi * t_grid * (1 - t_grid)
    })

  } else if (warping == "flexible") {
    # Create multiple warping functions with incrementing seeds
    grid_warped <- sapply(1:n, function(i) {
      # Calculate seed for each warping function if a base seed is provided
      warp_seed <- if(is.null(seed)) NULL else seed + i

      # Generate and apply warping function
      generate_warping_function(
        n_points = warping_points,
        gamma = warping_gamma,
        seed = warp_seed
      )(t_grid)
    })
  } else {
    stop("Warping type not specified correctly. Must be 'simple' or 'flexible'.")
  }

  # Transpose to get matrix with rows = curves, cols = points
  grid_warped <- t(grid_warped)

  # Evaluate base function at warped grid points and convert to functional data object
  warped_values <- tf::tfd(g(grid_warped), t_grid)

  # Calculate function variance for noise scaling
  fun_var <- tf::tf_fvar(warped_values, t_grid)

  # Create covariance matrix for noise (diagonal matrix with variance proportional to function variance)
  cov_mat <- matrix(0, n, n)
  diag(cov_mat) <- fun_var * noise_to_signal

  # Generate multivariate normal noise
  mu <- rep(0, n)
  noise <- MASS::mvrnorm(num_points, mu = mu, Sigma = cov_mat)

  # Convert noise to functional data object
  noise_tf <- tf::tfd(t(noise), t_grid)

  # Combine amplitude variation, warped function values, and noise
  curves <- warped_values * a_params + noise_tf

  # Return results
  return(list(
    curves = curves,
    t_grid = t_grid,
    base_function = g,
    amplitude_params = a_params,
    grid_warped = grid_warped,
    noise_grid = noise_tf
  ))
}


