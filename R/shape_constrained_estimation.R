#' Shape Constrained Function Estimation
#'
#' @description Performs function estimation constrained by peak and valley locations
#'
#' @param curve_data Functional data object of class tf
#' @param peak_locs Numeric vector of peak locations
#' @param valley_locs Numeric vector of valley locations
#' @param significant_peaks Numeric vector of significant peak indices
#' @param peak_labels Numeric vector of peak labels
#' @param mean_function Functional data object representing the mean function
#' @param basis_dim Integer specifying the basis dimension. Default is 8.
#' @param t_grid Numeric vector representing the time grid
#' @param rho Numeric smoothness penalty parameter. Default is 1e-9.
#' @param basis_type Character specifying basis type (default: "cr" for cubic regression splines)
#'
#' @return A functional data object of class tf representing the estimated function
#'
#' @import mgcv
#' @importFrom pracma fmincon
#' @import tf
#' @export
shape_constrained_estimation <- function(curve_data, peak_locs, valley_locs,
                                         significant_peaks, peak_labels,
                                         mean_function,
                                         basis_dim = 8, t_grid, rho = 1e-9,
                                         basis_type = "cr") {
  # Input validation
  if (missing(curve_data)) {
    stop("'curve_data' is required")
  }
  if (missing(peak_locs)) {
    stop("'peak_locs' is required")
  }
  if (missing(valley_locs)) {
    stop("'valley_locs' is required")
  }
  if (missing(significant_peaks)) {
    stop("'significant_peaks' is required")
  }
  if (missing(peak_labels)) {
    stop("'peak_labels' is required")
  }
  if (missing(mean_function)) {
    stop("'mean_function' is required")
  }
  if (missing(t_grid)) {
    stop("'t_grid' is required")
  }

  # Type checking
  if (!inherits(curve_data, "tf")) {
    stop("'curve_data' must be a tf object")
  }
  if (!is.numeric(peak_locs)) {
    stop("'peak_locs' must be a numeric vector")
  }
  if (!is.numeric(valley_locs)) {
    stop("'valley_locs' must be a numeric vector")
  }
  if (!is.numeric(significant_peaks)) {
    stop("'significant_peaks' must be a numeric vector")
  }
  if (!is.numeric(peak_labels)) {
    stop("'peak_labels' must be a numeric vector")
  }
  if (!inherits(mean_function, "tf")) {
    stop("'mean_function' must be a tf object")
  }
  if (!is.numeric(basis_dim) || length(basis_dim) != 1 || basis_dim != round(basis_dim) || basis_dim <= 0) {
    stop("'basis_dim' must be a positive integer")
  }
  if (!is.numeric(t_grid)) {
    stop("'t_grid' must be a numeric vector")
  }
  if (!is.numeric(rho) || length(rho) != 1 || rho < 0) {
    stop("'rho' must be a non-negative numeric value")
  }
  if (!is.character(basis_type) || length(basis_type) != 1) {
    stop("'basis_type' must be a character string")
  }

  # Check for required packages
  if (!requireNamespace("mgcv", quietly = TRUE)) {
    stop("Package 'mgcv' is needed for this function to work. Please install it.")
  }
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("Package 'pracma' is needed for this function to work. Please install it.")
  }
  if (!requireNamespace("tf", quietly = TRUE)) {
    stop("Package 'tf' is needed for this function to work. Please install it.")
  }

  # Additional checks
  if (length(peak_locs) == 0) {
    stop("'peak_locs' cannot be empty")
  }
  if (length(peak_labels) != length(peak_locs)) {
    stop("'peak_labels' must have the same length as 'peak_locs'")
  }
  if (length(significant_peaks) == 0) {
    stop("'significant_peaks' cannot be empty")
  }
  if (length(t_grid) < 2) {
    stop("'t_grid' must have at least two points")
  }

  # Check for reconstruct_function
  if (!exists("reconstruct_function")) {
    stop("Function 'reconstruct_function' is required but not available")
  }

  # Filter for significant peaks
  names(peak_locs) <- peak_labels
  peak_locs <- peak_locs[as.character(significant_peaks)]

  # Check bounds
  if (mean_function[, t_grid[2]] > mean_function[, t_grid[2]]) {
    peak_locs <- c(peak_locs, t_grid[2])
  } else {
    valley_locs <- c(valley_locs, t_grid[2])
  }
  if (mean_function[, t_grid[length(t_grid) - 1]] >
      mean_function[, t_grid[length(t_grid) - 2]]) {
    peak_locs <- c(peak_locs, t_grid[length(t_grid) - 1])
  } else {
    valley_locs <- c(valley_locs, t_grid[length(t_grid) - 1])
  }

  # Find deepest valley between consecutive peaks
  peak_locs <- sort(peak_locs)
  valley_locs <- sort(valley_locs)
  split_points <- findInterval(valley_locs, peak_locs)

  tryCatch({
    valley_locs <- tapply(valley_locs, split_points, function(v) {
      v[which.min(mean_function[, v])]
    })
  }, error = function(e) {
    stop("Error finding deepest valleys: ", e$message)
  })

  idx <- sort(unique(c(peak_locs, valley_locs)))
  heights <- mean_function[, idx]
  n_extrema <- length(idx)

  # Create the smooth specification
  s_spec <- tryCatch({
    mgcv::s(t, k = basis_dim, bs = basis_type, m = 2)
  }, error = function(e) {
    stop("Error creating smooth specification: ", e$message)
  })

  # Create the basis constructor object with proper data structure
  construction_data <- data.frame(
    t = t_grid,
    .name_repair = "minimal"
  )

  # Construct the smooth object using mgcv's machinery
  smooth_object <- tryCatch({
    mgcv::smooth.construct(
      s_spec,
      data = construction_data,
      knots = NULL
    )
  }, error = function(e) {
    stop("Error constructing smooth object: ", e$message)
  })

  # Create inequality constraint matrix A and vector b
  # This ensures alternating maxima and minima
  A_ineq <- matrix(0, nrow = n_extrema-1, ncol = n_extrema)
  for (i in 1:(n_extrema-1)) {
    A_ineq[i, i] <- 1
    A_ineq[i, i+1] <- -1
  }

  # Determine if starting with maximum or minimum
  start <- peak_locs[1] == idx[1]

  # Adjust signs based on whether we start with maximum or minimum
  if (start) {
    # Select odd-numbered rows (1, 3, 5, ...) and multiply them by -1
    A_ineq[seq(1, nrow(A_ineq), by = 2), ] <- -A_ineq[seq(1, nrow(A_ineq), by = 2), ]
  } else {
    # Select even-numbered rows (2, 4, 6, ...) and multiply them by -1
    A_ineq[seq(2, nrow(A_ineq), by = 2), ] <- -A_ineq[seq(2, nrow(A_ineq), by = 2), ]
  }

  # Add zeros for basis coefficients to constraint matrix
  A_ineq <- cbind(matrix(0, nrow = n_extrema-1, ncol = basis_dim), A_ineq)
  b_ineq <- rep(-1e-6, n_extrema-1)

  # Define objective function
  obj_fn <- function(param) {
    # Reconstruct function based on parameters
    f_reconstr <- tryCatch({
      reconstruct_function(param = param, t_grid = t_grid,
                           idx = idx,
                           construct_data = construction_data,
                           smooth_object = smooth_object)
    }, error = function(e) {
      warning("Error in reconstruct_function: ", e$message)
      return(Inf)  # Return Inf to signal a bad parameter set
    })

    # Squared residuals
    residuals <- f_reconstr - curve_data

    sq_res_integral <- tryCatch({
      tf::tf_integrate(residuals^2, lower = tf::tf_domain(residuals)[1],
                       upper = tf::tf_domain(residuals)[2])
    }, error = function(e) {
      warning("Error integrating residuals: ", e$message)
      return(Inf)
    })

    # Squared second derivative as smoothness penalty
    second_deriv <- tryCatch({
      tf::tf_derive(f_reconstr, tf::tf_arg(f_reconstr), order = 2)
    }, error = function(e) {
      warning("Error calculating second derivative: ", e$message)
      return(Inf)
    })

    second_deriv_sq <- second_deriv^2

    smoothness_penalty <- tryCatch({
      tf::tf_integrate(second_deriv_sq,
                       lower = tf::tf_domain(second_deriv_sq)[1],
                       upper = tf::tf_domain(second_deriv_sq)[2])
    }, error = function(e) {
      warning("Error integrating second derivative: ", e$message)
      return(Inf)
    })

    # Objective: combination of fit and smoothness
    obj <- sqrt(sum(sq_res_integral, na.rm = TRUE)) +
      rho * smoothness_penalty

    return(obj)
  }

  # Initial values
  init_val <- c(rep(0, basis_dim), heights)

  # Run optimization using fmincon
  result <- tryCatch({
    pracma::fmincon(
      x0 = init_val,                # Initial values
      fn = obj_fn,                  # Objective function
      A = A_ineq,                   # Inequality constraint matrix
      b = b_ineq,                   # Inequality constraint vector
      method = "SQP",               # Sequential Quadratic Programming (like MATLAB)
      tol = 1e-6,                   # Tolerance for convergence
      maxfeval = 100000,             # Maximum function evaluations
      maxiter = 50000
    )
  }, error = function(e) {
    stop("Optimization failed: ", e$message)
  })

  # Reconstruct final function
  param_est <- result$par

  fn_est <- tryCatch({
    reconstruct_function(param_est, t_grid, idx, construction_data, smooth_object)
  }, error = function(e) {
    stop("Error reconstructing final function: ", e$message)
  })

  return(fn_est)
}


