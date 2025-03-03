#' Create Parameter Grid for Benchmarking
#'
#' @description
#' Creates a data frame containing all combinations of parameter settings for benchmarking.
#' Handles hierarchical parameters by allowing nested lists.
#'
#' @param ... Named parameters where each parameter can be a vector of values to test.
#'            For hierarchical parameters, use lists with nested structure.
#' @param random_subset Numeric, if provided, randomly selects this many rows from the full grid.
#' @param seed Integer, random seed for reproducibility when using random_subset.
#'
#' @return A data frame with each row representing a unique combination of parameters.
#'
#' @examples
#' # Simple grid
#' grid <- create_benchmark_grid(
#'   n = c(25, 50),
#'   lambda = c(0.1, 0.5, 1.0),
#'   warping = c("simple", "flexible")
#' )
#'
#' # With hierarchical parameters
#' grid <- create_benchmark_grid(
#'   n = c(25, 50),
#'   warping = list(
#'     simple = list(),
#'     flexible = list(
#'       warping_gamma = c(1, 5),
#'       warping_points = c(3, 5)
#'     )
#'   )
#' )
create_benchmark_grid <- function(..., random_subset = NULL, seed = NULL) {
  # Get the parameters and their values
  params <- list(...)

  # Check if any parameters are hierarchical (lists)
  hierarchical_params <- sapply(params, is.list)

  if (any(hierarchical_params)) {
    # Handle hierarchical parameters
    flat_grid <- expand_hierarchical_params(params)
  } else {
    # Simple case: use expand.grid directly
    flat_grid <- expand.grid(params, stringsAsFactors = FALSE)
  }

  # Take random subset if requested
  if (!is.null(random_subset)) {
    if (!is.null(seed)) set.seed(seed)
    if (random_subset < nrow(flat_grid)) {
      flat_grid <- flat_grid[sample(nrow(flat_grid), random_subset), ]
    }
  }

  # Add an ID column for tracking
  flat_grid$config_id <- 1:nrow(flat_grid)

  return(flat_grid)
}

#' Expand Hierarchical Parameters
#'
#' @param params A list of parameters where some may be hierarchical.
#'
#' @return A data frame with expanded parameter combinations.
#'
#' @keywords internal
expand_hierarchical_params <- function(params) {
  # Separate simple and hierarchical parameters
  simple_params <- params[!sapply(params, is.list)]
  hier_params <- params[sapply(params, is.list)]

  # Create base grid from simple parameters
  if (length(simple_params) > 0) {
    base_grid <- expand.grid(simple_params, stringsAsFactors = FALSE)
  } else {
    base_grid <- data.frame(dummy = 1)
  }

  # Process each hierarchical parameter
  result_grid <- base_grid

  for (param_name in names(hier_params)) {
    param_values <- hier_params[[param_name]]
    temp_grids <- list()

    # Collect all possible sub-parameter names across all values
    all_sub_params <- unique(unlist(lapply(param_values, function(x) names(x))))

    # For each value of the hierarchical parameter
    for (value_name in names(param_values)) {
      # Get the sub-parameters for this value
      sub_params <- param_values[[value_name]]

      if (length(sub_params) == 0) {
        # No sub-parameters, just add the value
        temp_grid <- result_grid
        temp_grid[[param_name]] <- value_name

        # Add all possible sub-parameters with NA values
        for (sub_param in all_sub_params) {
          if (!(sub_param %in% names(temp_grid))) {
            temp_grid[[sub_param]] <- NA
          }
        }

        temp_grids[[length(temp_grids) + 1]] <- temp_grid
      } else {
        # Expand the sub-parameters
        sub_grid <- expand.grid(sub_params, stringsAsFactors = FALSE)

        # Combine with the current result grid
        for (i in 1:nrow(sub_grid)) {
          temp_grid <- result_grid
          temp_grid[[param_name]] <- value_name

          # Add sub-parameter columns
          for (sub_param in names(sub_params)) {
            temp_grid[[sub_param]] <- sub_grid[i, sub_param]
          }

          # Add missing sub-parameters with NA values
          for (sub_param in setdiff(all_sub_params, names(sub_params))) {
            if (!(sub_param %in% names(temp_grid))) {
              temp_grid[[sub_param]] <- NA
            }
          }

          temp_grids[[length(temp_grids) + 1]] <- temp_grid
        }
      }
    }

    # Ensure all grids have the same columns before combining
    all_cols <- unique(unlist(lapply(temp_grids, names)))

    # Add any missing columns to each grid with NA values
    for (i in seq_along(temp_grids)) {
      missing_cols <- setdiff(all_cols, names(temp_grids[[i]]))
      for (col in missing_cols) {
        temp_grids[[i]][[col]] <- NA
      }
    }

    # Combine all temporary grids
    result_grid <- do.call(rbind, temp_grids)
  }

  # Remove dummy column if it exists
  if ("dummy" %in% names(result_grid)) {
    result_grid$dummy <- NULL
  }

  return(result_grid)
}
