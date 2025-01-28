library(pracma)


# Function for shape constrained estimation
shape_constrained_estimation <- function(X, g_template, peak_locs, valley_locs,
                                         K, t_grid, rho, basis_type = "Fourier") {
  # Find number of local maxima in template
  n_maxima <- length(peak_locs)

  if (n_maxima == 0) {
    return(list(fn_est = NULL, K = NULL))
  }

  # Create basis matrix
  b <- generate_basis(K, t_grid, basis_type)
  b_size <- ncol(b)

  # Find extrema
  # Code from paper does it in function, unclear if necessary
  # idxM <- find_local_maxima(g_template)
  # idxm <- find_local_minima(g_template)

  temp_idx <- sort(c(peak_locs, valley_locs))
  heights <- g_template[, temp_idx]

  # Determine if starting with maximum or minimum
  start <- peak_locs[[1]] == temp_idx[[1]]

  n <- length(temp_idx)

  # Create inequality constraint matrix A and vector b
  # This ensures alternating maxima and minima
  A_ineq <- matrix(0, nrow = n-1, ncol = n)
  for (i in 1:(n-1)) {
    A_ineq[i, i] <- 1
    A_ineq[i, i+1] <- -1
  }

  # Adjust signs based on whether we start with maximum or minimum
  if (start) {
    A_ineq[, seq(1, n, by = 2)] <- -A_ineq[, seq(1, n, by = 2)]
  } else {
    A_ineq[, seq(2, n, by = 2)] <- -A_ineq[, seq(2, n, by = 2)]
  }

  # Add zeros for basis coefficients to constraint matrix
  A_ineq <- cbind(matrix(0, nrow = n-1, ncol = b_size), A_ineq)
  b_ineq <- rep(-1e-6, n-1)

  # Objective function
  obj_fn <- function(param) {
    beta <- param[1:b_size]
    h <- param[(b_size+1):length(param)]

    # Reconstruct function
    f_reconstr <- b %*% beta
    for (i in seq_along(temp_idx)) {
      f_reconstr[temp_idx[i]] <- h[i]
    }

    # Calculate objective
    #TODO: keep TF?
    #TODO:
    X <- X[, t_grid]
    residuals <- X - matrix(rep(f_reconstr, each = nrow(X)), nrow = nrow(X))

    # Calculate second derivative
    first_grad <- diff(f_reconstr) / diff(t_grid)
    second_grad <- diff(first_grad) / diff(t_grid[-1])
    obj <- sqrt(sum(residuals^2, na.rm = TRUE)) + rho * sum(second_grad^2, na.rm = TRUE)
    return(obj)
  }

  # Initial values
  init_val <- c(rep(1e-10, b_size), heights)

  # Run optimization using fmincon
  result <- fmincon(
    x0 = init_val,                # Initial values
    fn = obj_fn,                  # Objective function
    A = A_ineq,                   # Inequality constraint matrix
    b = b_ineq,                   # Inequality constraint vector
    method = "SQP",               # Sequential Quadratic Programming (like MATLAB)
    tol = 1e-6,                   # Tolerance for convergence
    maxfeval = 20000,             # Maximum function evaluations
    maxiter = 10000               # Maximum iterations
  )

  # Reconstruct final function
  param_est <- result$par
  beta_est <- param_est[1:b_size]
  h_est <- param_est[(b_size+1):length(param_est)]

  fn_est <- b %*% beta_est
  for (i in seq_along(temp_idx)) {
    fn_est[temp_idx[i]] <- h_est[i]
  }

  fn_est <- as.vector(fn_est)
  fn_tf <- tfd(fn_est, t_grid)

  return(fn_tf)
}
