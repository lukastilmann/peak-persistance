#' Align Functional Data to a Target Function
#'
#' @description
#' Aligns a set of functional data objects to a target function (typically a mean function)
#' using the Fisher-Rao metric and square-root velocity framework. This function serves
#' as a wrapper for the `fdasrvf::multiple_align_functions` function with additional
#' input validation and convenient return format.
#'
#' @param function_curves A functional data object inheriting from class `tf` containing
#'        the curves to be aligned.
#' @param target_function A functional data object inheriting from class `tf` representing
#'        the target function (typically a mean function) to align to. Must contain exactly
#'        one function.
#' @param t_grid Numeric vector of time points for evaluation. If NULL (default),
#'        the function will use the evaluation grid from the target_function.
#' @param lambda Numeric value specifying the regularization parameter.
#'        Controls the smoothness of the warping functions.
#'
#' @return A list containing two elements:
#'   \item{aligned_functions}{A `tfd` object containing the aligned curves}
#'   \item{warping_functions}{A `tfd` object containing the warping functions}
#'
#' @details
#' This function aligns a set of functions to a target function (typically a mean or template)
#' using the Fisher-Rao metric. The alignment is performed through time warping
#' functions that optimize the similarity between each curve and the target while
#' maintaining smoothness controlled by the `lambda` parameter.
#'
#' Unlike the general alignment process where the mean is iteratively updated,
#' this function aligns each curve to a fixed target function.
#'
#' @examples
#' \dontrun{
#' # Create some sample functional data
#' data <- tfd(matrix(rnorm(500), nrow=5), arg=seq(0, 1, length.out=100))
#'
#' # Calculate the mean function
#' mean_func <- tf_mean(data)
#'
#' # Align all curves to the mean
#' aligned <- align_to_mean(data, mean_func, lambda=0.1)
#'
#' # Plot original vs aligned curves
#' par(mfrow=c(1,2))
#' plot(data, main="Original curves")
#' plot(aligned$aligned_functions, main="Aligned curves")
#' }
#'
#' @export
align_to_mean <- function(function_curves, target_function, t_grid = NULL, lambda = 0.1) {
  # Input validation for function_curves
  if (!inherits(function_curves, "tf")) {
    stop("function_curves must inherit from class 'tf'")
  }

  # Input validation for target_function
  if (!inherits(target_function, "tf")) {
    stop("target_function must inherit from class 'tf'")
  }

  # Check that target_function contains exactly one function
  if (length(target_function) != 1) {
    stop("target_function must contain exactly one function")
  }

  # If t_grid is NULL, get it from target_function
  if (is.null(t_grid)) {
    t_grid <- tf_arg(target_function)
    # Since target_function is a single function, tf_arg should return a simple vector
    # If it's still a list for some reason, provide a clear error
    if (is.list(t_grid)) {
      stop("Unexpected list returned from tf_arg(target_function). The target function may have an irregular grid structure.")
    }
  } else {
    # Validate t_grid if provided
    if (!is.numeric(t_grid) || any(is.na(t_grid)) || any(diff(t_grid) <= 0)) {
      stop("t_grid must be a numeric vector with strictly increasing values")
    }
  }

  # Validate lambda
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("lambda must be a non-negative numeric value")
  }

  # Getting data into right format
  functions <- t(function_curves[, t_grid])
  mu <- as.vector(target_function[, t_grid])

  # Align to target function
  res <- fdasrvf::multiple_align_functions(functions, t_grid, mu,
                                           showplot = FALSE, lambda = lambda)

  # Aligned functions and warping functions in tf format
  aligned_functions <- t(res$fn)
  aligned_functions_tf <- tfd(aligned_functions, t_grid)
  warping_functions <- t(res$warping_functions)
  warping_functions_tf <- tfd(warping_functions, t_grid)

  return(list(
    aligned_functions = aligned_functions_tf,
    warping_functions = warping_functions_tf
  ))
}
