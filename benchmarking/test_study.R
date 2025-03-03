devtools::load_all(".") # do not "source" package code - load the package!
library(tidyfun) # so geom_spaghetti is available
source("./benchmarking/create_settings.R")
source("./benchmarking/run_simulation.R")
source("./benchmarking/calculate_benchmark_metrics.R")
source("./benchmarking/benchmark_study.R")


fun_1 <- generate_benchmark_function(list("normal", "normal"),
                                     c(10, 10), c(0.3, 0.7), c(10, 10))

fun_2 <- generate_benchmark_function(
  bases = list("normal", "normal", "normal", "normal"),  # x^1, x^2, normal, sine
  coefficients = list(
    10,
    5,
    3,
    7# Coefficient for sine
  ),
  shifts = c(0.2, 0.5, 0.7, 0.9),     # Shifts for each basis
  scale_factors = c(15, 10, 20, 10)      # Amplitude factors
)

functions_list <- c(fun_1, fun_2)

# With hierarchical parameters
grid <- create_benchmark_grid(
  g_id = c(1, 2),
  n_lambda = c(25),
  n = c(20, 50),
  lambda_search_threshold = c(0.15),
  n_points = c(100, 20),
  warping = list(
    simple = list(),
    flexible = list(
      warping_gamma = c(5),
      warping_points = c(5)
    )
  ),
  lambda_search_min_bound = c(0.1),
  penalty = c("roughness", "geodesic"),
  pers_method = c("clustering", "threshold")
)


benchmark_study <- run_benchmark_study(grid, functions_list,
                                       output_dir = "./benchmarking/results/test",
                                       save_plots = TRUE, runs_per_config = 3)
