devtools::load_all(".") # do not "source" package code - load the package!
library(tidyfun) # so geom_spaghetti is available
source("./benchmarking/create_settings.R")
source("./benchmarking/run_simulation.R")
source("./benchmarking/calculate_benchmark_metrics.R")
source("./benchmarking/benchmark_study.R")


fun_1 <- generate_benchmark_function(list("normal", "normal"),
                                     c(10, 10), c(0.3, 0.7), c(10, 10))

fun_2 <- generate_benchmark_function(
  bases = list("normal", "normal", "normal"),
  coefficients = list(
    10,
    10,
    10
  ),
  shifts = c(0.25, 0.5, 0.75),     # Shifts for each basis
  scale_factors = c(15, 15, 15)      # Amplitude factors
)

functions_list <- c(fun_1, fun_2)

# With hierarchical parameters
grid <- create_benchmark_grid(
  g_id = c(1, 2),
  n_lambda = c(5),
  lambda_search_min_bound = c(0.5),
  warping = list(
    simple = list(),
    flexible = list(
      warping_gamma = c(1, 5),
      warping_points = c(1, 5)
    )
  )
)

grid <- grid[1:5,]

benchmark_study <- run_benchmark_study(grid, functions_list,
                                       output_dir = "./benchmarking/results/test_3",
                                       save_plots = TRUE, seed = 123)
