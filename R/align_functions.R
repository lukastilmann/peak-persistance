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

    # Get the initial values of original functions
    f0_init <- data_aligned_tw$f0[1,]
    #View(data_aligned_tw)

    # Calculate the mean/median that was used
    if(data_aligned_tw$call$centroid_type == "mean") {
      center_value <- mean(f0_init)
    } else {
      center_value <- median(f0_init)
    }

    # Get the aligned functions
    fn <- data_aligned_tw$fn

    # Recover each original function
    recovered <- matrix(0, nrow=nrow(fn), ncol=ncol(fn))
    for(i in 1:ncol(fn)) {
      # Add back the difference in initial values
      offset <- f0_init[i] - center_value
      recovered[,i] <- fn[,i] + offset
    }

    result <- t(recovered)
  } else {
    print("csa for curve_srvf_align or tw for time_warping implemented")
  }

  aligned_curves <- tfd(result, arg = t_grid)
  return(aligned_curves)
}
