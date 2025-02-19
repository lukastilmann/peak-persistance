library(fdasrvf)

align_functions <- function(fun_curves, lambda = 0.0, parallel = FALSE,
                            max_iter = 10,
                            t_grid = seq(0, 1, length.out = 100),
                            verbosity = "low"){

  curves <- fun_curves$curves
  curves_arr_tw <- t(curves[, t_grid])

  # Flag that notes convergence (not perfectly, since time_warping doesn't mark
  # convergence, but just throws a warning when last iteration is entered, where
  # it may converge)
  converged <- TRUE

  # Clean up any lingering parallel connections
  closeAllConnections()

  # If that's not enough, try this more aggressive cleanup:
  if (exists("cl")) {
    try(parallel::stopCluster(cl), silent = TRUE)
  }

  alignment_result <- withCallingHandlers(
    {
      message(paste("Alignment using lambda =", lambda))
      time_warping(curves_arr_tw, t_grid, lambda = lambda, parallel = parallel,
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

  return(list(
    aligned_curves = aligned_curves,
    converged = converged
  ))
}

