#' Reconstruct Function from Parameters
#'
#' @description
#' Reconstructs a functional data object from optimization parameters, combining
#' basis coefficients and height values at extrema.
#'
#' @param param Numeric vector of parameters (basis coefficients and heights)
#' @param t_grid Numeric vector of time grid points
#' @param idx Numeric vector of indices for extrema points
#' @param construct_data Data frame containing time grid data
#' @param smooth_object Smooth object from mgcv
#'
#' @return A functional data object of class tf
#'
#' @importFrom mgcv PredictMat
#' @importFrom stats splinefun
#' @import tf
reconstruct_function <- function(param, t_grid, idx, construct_data, smooth_object) {
  # Input validation
  if (missing(param)) {
    stop("'param' is required")
  }
  if (missing(t_grid)) {
    stop("'t_grid' is required")
  }
  if (missing(idx)) {
    stop("'idx' is required")
  }
  if (missing(construct_data)) {
    stop("'construct_data' is required")
  }
  if (missing(smooth_object)) {
    stop("'smooth_object' is required")
  }

  # Type checking
  if (!is.numeric(param)) {
    stop("'param' must be a numeric vector")
  }
  if (!is.numeric(t_grid)) {
    stop("'t_grid' must be a numeric vector")
  }
  if (!is.numeric(idx)) {
    stop("'idx' must be a numeric vector")
  }
  if (!is.data.frame(construct_data)) {
    stop("'construct_data' must be a data frame")
  }
  #TODO: input check for smooth_object

  # Check parameters length
  basis_size <- ncol(smooth_object$X)
  if (length(param) <= basis_size) {
    stop("'param' is too short: must contain basis coefficients and heights")
  }

  # Check idx elements are within t_grid range
  if (any(idx < min(t_grid) | idx > max(t_grid))) {
    stop("All elements in 'idx' must be within the range of 't_grid'")
  }

  # Reconstruct function from parameters
  beta <- param[1:basis_size]
  heights <- param[(basis_size+1):length(param)]

  # Check heights length matches idx length
  if (length(heights) != length(idx)) {
    stop("Number of heights does not match number of extrema points")
  }

  # Computing warping function from spline weights
  spline_fun <- tryCatch({
    mgcv::PredictMat(smooth_object, construct_data) %*% beta
  }, error = function(e) {
    stop("Error computing spline function: ", e$message)
  })

  p <- rep(1, length(spline_fun))
  nv <- norm(spline_fun, type = "2") / length(spline_fun)

  # Avoid division by zero
  if (nv < 1e-10) {
    exponential_mapping <- p
  } else {
    exponential_mapping <- cos(nv) * p + sin(nv) * spline_fun / nv
  }

  # Warping time grid
  warping_fun <- cumsum(exponential_mapping^2) * (length(t_grid)-1)^-1
  t_grid_warped <- warping_fun / warping_fun[length(warping_fun)]

  # Create function using spline interpolation
  fun_est <- tryCatch({
    stats::splinefun(idx, heights, "monoH.FC")(t_grid_warped)
  }, error = function(e) {
    stop("Error in spline interpolation: ", e$message)
  })

  # Convert to tf object
  fun_est_tf <- tryCatch({
    tf::tfd(fun_est, t_grid)
  }, error = function(e) {
    stop("Error converting to tf object: ", e$message)
  })

  return(fun_est_tf)
}
