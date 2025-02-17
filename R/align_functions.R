library(fdasrvf)

align_functions <- function(fun_curves, lambda = 0.0, parallel = FALSE,
                            max_iter = 10,
                            t_grid = seq(0, 1, length.out = 100),
                            func = "tw"){

  if (func == "csa"){
    # data format needed for fdasrvf library
    curves_arr <- tf_to_fdasrvf(fun_curves, "curves", t_grid)

    # TODO: test which arguments lead to crashes and whether parallelizing
    # here is efficient
    data_aligned <- curve_srvf_align(curves_arr, mode = "O", rotated = TRUE,
                                     maxit = max_iter, parallel = parallel,
                                     lambda = lambda,
                                     scale = FALSE)

    # Then convert back to tf format
    result <- fdasrvf_to_matrix(data_aligned, t_grid)
  } else if (func == "tw"){
    # Using the time_warping() function
    curves <- fun_curves$curves
    curves_arr_tw <- t(curves[, t_grid])
    # TODO: any way to turn down verbosity?
    data_aligned_tw <- time_warping(curves_arr_tw, t_grid, lambda = lambda,
                                    parallel = TRUE, max_iter = 20)

    # Get the aligned functions
    fn <- data_aligned_tw$fn
    result <- t(fn)
  } else {
    print("csa for curve_srvf_align or tw for time_warping implemented")
  }

  aligned_curves <- tfd(result, arg = t_grid)
  return(aligned_curves)
}

