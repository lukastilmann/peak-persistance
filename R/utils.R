
tf_to_fdasrvf <- function(tf_data, tf_col, grid = NULL) {
  # Convert tf object to wide format using tf_spread with specified grid
  #wide_data <- tf_spread(tf_data, !!enquo(tf_col), arg = grid)
  wide_data <- tf_spread(tf_data, tf_col, arg = grid)

  # Get the column names that contain the function evaluations
  #value_cols <- grep(paste0(quo_name(enquo(tf_col)), "_"), names(wide_data), value = TRUE)
  value_cols <- grep(paste0(tf_col , "_"), names(wide_data), value = TRUE)

  # Extract just the function values into a matrix
  func_matrix <- as.matrix(wide_data[, value_cols])

  # Get dimensions
  N <- nrow(func_matrix)  # number of curves
  T <- ncol(func_matrix)  # number of time points
  n <- 2  # dimension (R²)

  # Create the 3D array with proper dimensions: n × T × N
  beta <- array(0, dim = c(n, T, N))

  # Fill the array:
  # First coordinate (x) is the grid points
  # Second coordinate (y) is the function values
  beta[1, , ] <- matrix(grid, nrow = 1, ncol = T) %>%
    replicate(n = N, simplify = "array")
  beta[2, , ] <- t(func_matrix)

  return(beta)
}


fdasrvf_to_matrix <- function(aligned_curves, grid) {
  # Extract aligned curves from fdasrvf output (betan)
  beta_aligned <- aligned_curves$betan

  # Extract y-coordinates (dimension 2)
  f <- beta_aligned[2, , ]

  #TODO: fdasrvf seems to standardize to 0 mean
  # figure out how to add back the mean
  # centers <- aligned_curves$cent
  #
  # for(i in 1:ncol(f)) {
  #   f[, i] <- f[, i] + centers[2, i]
  # }

  # Return M x N matrix (transpose so time points are rows)
  t(f)
}


get_lambda <- function(col_name) {
  as.numeric(sub("aligned_lambda_", "", col_name))
}
