#' Align a Function to Another Using SRVF Methodology
#'
#' @description
#' Aligns one functional object to another using the Square Root Velocity Function (SRVF) framework.
#' This function performs temporal alignment (warping) to find the optimal time transformation
#' that makes the two functions as similar as possible while preserving shape characteristics.
#'
#' @param function_1 An object of class 'tf' representing the first function to be used as reference.
#'   Must contain exactly one function.
#' @param function_2 An object of class 'tf' representing the second function to be aligned to the first.
#'   Must contain exactly one function.
#' @param t_grid A numeric vector specifying the time grid over which to perform the alignment.
#'   If NULL (default), the time grid from function_1 will be used.
#' @param lambda A non-negative numeric value controlling the amount of regularization in the
#'   alignment process. Higher values result in smoother warping functions. Default is 0.
#'
#' @return A list containing four components:
#'   \item{function_2_aligned}{A 'tf' object representing function_2 after alignment to function_1}
#'   \item{warping_function}{A 'tf' object representing the warping function that transforms function_2}
#'   \item{q_1}{A 'tf' object containing the SRVF of function_1}
#'   \item{q_2_aligned}{A 'tf' object containing the SRVF of function_2 after alignment}
#'
#' @details
#' The function uses the Square Root Velocity Function (SRVF) framework from the fdasrvf package
#' to find optimal time-warping. The alignment process finds a warping function that minimizes
#' the distance between the SRVFs of the two functions, with an optional regularization term
#' controlled by lambda. The Dynamic Programming algorithm is used for optimization.
#'
#' @note
#' There appears to be a possible bug in the return statement, where 'q_2_tf' is referenced
#' but not defined. This might need to be corrected to 'q_2_aligned_tf'.
#'
#' @examples
#' \dontrun{
#' # Create two example tf objects
#' t <- seq(0, 1, length.out = 100)
#' f1 <- tfd(sin(2*pi*t), t)
#' f2 <- tfd(sin(2*pi*t + 0.5), t)
#'
#' # Align f2 to f1
#' result <- align_to_mean(f1, f2, lambda = 0.1)
#'
#' # Plot the results
#' plot(f1, col = "black", main = "Alignment Result")
#' lines(f2, col = "red")
#' lines(result$function_2_aligned, col = "blue")
#' }
#'
#' @seealso
#' \code{\link[fdasrvf]{f_to_srvf}}, \code{\link[fdasrvf]{warp_f_gamma}}
#'
#' @export
align_pairwise <- function(function_1, function_2, t_grid = NULL, lambda = 0) {
  # Input validation for function_curves
  if (!inherits(function_1, "tf")) {
    stop("function_curves must inherit from class 'tf'")
  }

  # Input validation for target_function
  if (!inherits(function_2, "tf")) {
    stop("target_function must inherit from class 'tf'")
  }

  # Check that target_function contains exactly one function
  if (length(function_1) != 1 || length(function_2) != 1) {
    stop("Both functions must be tf of length 1.")
  }

  # If t_grid is NULL, get it from function_1
  if (is.null(t_grid)) {
    t_grid <- tf_arg(function_1)
  }

  # Validate lambda
  if (!is.numeric(lambda) || length(lambda) != 1 || lambda < 0) {
    stop("lambda must be a non-negative numeric value")
  }

  # Getting data into right format
  fun_1 <- as.vector(function_1[, t_grid])
  fun_2 <- as.vector(function_2[, t_grid])

  # SRVF of both functions
  q_1 = fdasrvf::f_to_srvf(fun_1, t_grid)
  q_2 = fdasrvf::f_to_srvf(fun_2, t_grid)

  # Align to target function
  warping_function <- optimum.reparam(q_1, t_grid, q_2, t_grid, lambda, "roughness",
                         "DP", f1o=fun_1[1], f2o=fun_2[1])

  fun_2_aligned <- fdasrvf::warp_f_gamma(fun_2, t_grid, warping_function)
  q_2_aligned <- fdasrvf::warp_q_gamma(q_2, t_grid, warping_function)

  # Aligned function 2 and warping function as tfd
  function_2_aligned_tf <- tfd(fun_2_aligned, t_grid)
  warping_function_tf <- tfd(warping_function, t_grid)
  q_1_tf <- tfd(q_1, t_grid)
  q_2_aligned_tf <- tfd(q_2_aligned, t_grid)

  return(list(
    function_2_aligned = function_2_aligned_tf,
    warping_function = warping_function_tf,
    q_1 = q_1_tf,
    q_2_aligned = q_2_aligned_tf
  ))
}
