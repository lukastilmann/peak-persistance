#' Generate basis functions for functional data analysis
#'
#' This function creates different types of basis functions commonly used in
#' functional data analysis, particularly for shape-constrained estimation.
#'
#' @param K Parameter for basis functions:
#'          - For Fourier/Cosine: number of basis functions (scalar)
#'          - For Meyer: vector of length 2 specifying range
#' @param t Vector of time points or domain values
#' @param type Type of basis functions to generate: "Fourier", "Cosine", or "Meyer"
#' @return Matrix where each column is a basis function evaluated at points t
#'
#' @examples
#' t <- seq(0, 1, length.out = 100)
#' # Generate Fourier basis with 5 components
#' fourier_basis <- generate_basis(5, t, "Fourier")
#'
generate_basis <- function(K, t, type = "Fourier") {
  # Input validation
  if (!is.numeric(t)) {
    stop("Time points 't' must be numeric")
  }

  # Validate K based on basis type
  if (type %in% c("Fourier", "Cosine")) {
    if (length(K) != 1 || !is.numeric(K)) {
      stop("For Fourier or Cosine basis, K must be a scalar")
    }
  } else if (type == "Meyer") {
    if (length(K) != 2 || !is.numeric(K)) {
      stop("For Meyer basis, K must be a numeric vector of length 2")
    }
  } else {
    stop("Unknown basis type. Must be one of: 'Fourier', 'Cosine', or 'Meyer'")
  }

  # Generate basis functions based on type
  if (type == "Fourier") {
    # Create sequence of k values
    k_seq <- 1:K

    # Initialize basis matrix
    # We need 2K functions: K sine and K cosine terms
    b <- matrix(0, nrow = length(t), ncol = 2 * K)

    # Generate sine and cosine terms
    for (k in k_seq) {
      # Cosine terms (odd columns)
      b[, 2*k-1] <- sqrt(2) * cos(2 * pi * k * t)
      # Sine terms (even columns)
      b[, 2*k] <- sqrt(2) * sin(2 * pi * k * t)
    }

  } else if (type == "Cosine") {
    # Create sequence of k values
    k_seq <- 1:K

    # Initialize basis matrix
    b <- matrix(0, nrow = length(t), ncol = K)

    # Generate cosine terms
    for (k in k_seq) {
      b[, k] <- sqrt(2) * cos(pi * k * t)
    }

  } else if (type == "Meyer") {
    # Note: Meyer wavelet basis implementation would go here
    # This is a complex wavelet basis that would require additional functions
    stop("Meyer basis implementation not currently available in R version")
  }

  return(b)
}

#' Helper function to visualize basis functions
#'
#' @param basis_matrix Matrix of basis functions from generate_basis()
#' @param t Time points used to generate the basis
#' @param type Type of basis that was generated
#'
plot_basis <- function(basis_matrix, t, type) {
  # Convert to data frame for ggplot
  df <- data.frame(
    t = rep(t, ncol(basis_matrix)),
    value = as.vector(basis_matrix),
    basis = factor(rep(1:ncol(basis_matrix), each = length(t)))
  )

  # Create plot using base R
  plot(NULL, xlim = range(t), ylim = range(basis_matrix),
       xlab = "t", ylab = "Basis Function Value",
       main = paste(type, "Basis Functions"))

  for (i in 1:ncol(basis_matrix)) {
    lines(t, basis_matrix[, i], col = i)
  }

  legend("topright", legend = paste("Basis", 1:ncol(basis_matrix)),
         col = 1:ncol(basis_matrix), lty = 1)
}
