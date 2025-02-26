#' Align Functional Data Objects
#'
#' @description
#' Performs alignment of functional data by finding time warping functions
#' that minimize a fitting criterion. This function serves as a wrapper for
#' the `time_warping` function with additional input validation and
#' user-friendly messaging.
#'
#' @param fun_curves A functional data object (inheriting from class `tf`) or
#'        a data frame with exactly one column that inherits from `tf`.
#' @param lambda Numeric value specifying the regularization parameter.
#'        Controls the smoothness of the warping functions. Default is 0.0.
#' @param parallel Logical indicating whether to use parallel processing.
#'        Default is FALSE.
#' @param max_iter Integer specifying the maximum number of iterations.
#'        Default is 10.
#' @param penalty Character string specify which penalty is used for function
#'        alignment.
#' @param t_grid Numeric vector of time points for evaluation. If NULL (default),
#'        the function will use the evaluation grid from the input object.
#' @param verbosity Character string specifying the level of progress messages.
#'        Must be either "low" or "high". Default is "low".
#'
#' @return A list containing three elements:
#'   \item{aligned_curves}{A `tfd` object containing the aligned curves}
#'   \item{warping_functions}{A `tfd` object containing the warping functions}
#'   \item{converged}{Logical indicating whether the algorithm converged}
#'
#' @details
#' The function aligns functional data using the Fisher-Rao metric and
#' square-root velocity framework. The alignment is performed by finding
#' time warping functions that minimize a criterion balancing data fit and
#' smoothness of the warping functions. The parameter `lambda` controls
#' this trade-off.
#'
#' @examples
#' \dontrun{
#' # Create some sample functional data
#' data <- tfd(matrix(rnorm(500), nrow=5), arg=seq(0, 1, length.out=100))
#'
#' # Align the curves
#' aligned <- align_functions(data, lambda=0.1)
#'
#' # Plot original vs aligned curves
#' par(mfrow=c(1,2))
#' plot(data, main="Original curves")
#' plot(aligned$aligned_curves, main="Aligned curves")
#' }
#'
#' @export
align_functions <- function(fun_curves, lambda = 0.0, parallel = FALSE,
                            max_iter = 10,
                            penalty = "roughness",
                            t_grid = NULL,
                            verbosity = "low") {
  # Input validation
  # Check if fun_curves is a tf object or a data frame with a tf column
  if (inherits(fun_curves, "tf")) {
    curves <- fun_curves
  } else if (is.data.frame(fun_curves)) {
    # Check if data frame has exactly one column
    if (ncol(fun_curves) != 1) {
      stop("If fun_curves is a data frame, it must have exactly one column.")
    }
    # Check if that column is a tf object
    if (!inherits(fun_curves[[1]], "tf")) {
      stop("The column in fun_curves must inherit from class 'tf'.")
    }
    curves <- fun_curves[[1]]
  } else {
    stop("fun_curves must be a 'tf' object or a data frame with one 'tf' column.")
  }

  # Validate verbosity
  if (!verbosity %in% c("low", "high")) {
    stop("verbosity must be either 'low' or 'high'.")
  }

  # Validate lambda
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("lambda must be a non-negative numeric value.")
  }

  # Validate max_iter
  if (!is.numeric(max_iter) || length(max_iter) != 1 ||
      max_iter != round(max_iter) || max_iter <= 0) {
    stop("max_iter must be a positive integer.")
  }

  # Validate penalty
  if (!penalty %in% c("roughness", "geodesic", "norm")) {
    stop("penalty must be either 'roughness', 'geodesic', or 'norm'.")
  }

  # Validate parallel
  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("parallel must be a logical value (TRUE or FALSE).")
  }

  # Use the evaluation grid from the input object if t_grid is NULL
  if (is.null(t_grid)) {
    t_grid <- tf_arg(curves)
    if (is.list(t_grid)) {
      # If it's a list (which could happen for irregular data),
      # use the first element's grid
      t_grid <- t_grid[[1]]
    }
  } else if (!is.numeric(t_grid) || any(is.na(t_grid)) ||
             any(diff(t_grid) <= 0)) {
    stop("t_grid must be a numeric vector with strictly increasing values.")
  }

  # Convert curves to array format required by time_warping
  curves_arr_tw <- t(curves[, t_grid])

  # Flag that notes convergence
  converged <- TRUE

  # Clean up any lingering parallel connections
  #closeAllConnections()

  alignment_result <- withCallingHandlers(
    {
      message(paste("Alignment using lambda =", lambda))
      fdasrvf::time_warping(curves_arr_tw, t_grid, lambda = lambda, parallel = parallel,
                   max_iter = max_iter)
    },
    message = function(m) {
      if (grepl("Computing Karcher .* of .* functions in SRSF space", m$message) ||
          grepl("Initializing...", m$message) ||
          grepl("Using lambda = ", m$message)) {
        # Suppress these messages completely
        invokeRestart("muffleMessage")
      }
      else if (grepl("maximal number of iterations", m$message)) {
        converged <<- FALSE
        invokeRestart("muffleMessage")
      }
      else if (grepl("Entering iteration ([0-9]+)", m$message)) {
        # Extract iteration number
        iter_num <- as.numeric(gsub(".*iteration ([0-9]+).*", "\\1", m$message))
        # Handle different verbosity levels
        if (verbosity == "low" && iter_num %% 2 == 0) {
          message(sprintf("Iteration: %d", iter_num))
          invokeRestart("muffleMessage")
        }
        else if (verbosity == "high") {
          message(sprintf("Iteration: %d", iter_num))
          invokeRestart("muffleMessage")
        }
        else if (verbosity == "low") {
          # Suppress odd-numbered iterations in low verbosity
          invokeRestart("muffleMessage")
        }
      }
    }
  )

  # Get the aligned functions
  fn <- alignment_result$fn
  aligned_curves <- tfd(t(fn), arg = t_grid)
  warp_fn <- alignment_result$warping_functions
  warping_functions <- tfd(t(warp_fn), arg = t_grid)

  return(list(
    aligned_curves = aligned_curves,
    warping_functions = warping_functions,
    converged = converged
  ))
}

