devtools::load_all(".")
pak::pak("tidyfun/tidyfun")
# Before setting up the cluster
devtools::install(quick = TRUE, quiet = TRUE)
library(tidyfun)
library(ggplot2)
library(foreach)
library(doParallel)
source("./benchmarking/create_settings.R")
source("./benchmarking/run_simulation.R")
source("./benchmarking/calculate_benchmark_metrics.R")
source("./benchmarking/benchmark_study.R")
source("./benchmarking/benchmark_utils.R")

t_grid <- seq(0, 1, length.out = 100)
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


fun_3 <- generate_benchmark_function(
  bases = list("sine", "normal", "normal", "normal", "normal"),
  coefficients = list(
    10,
    -30,
    -30,
    -30,
    -30
  ),
  shifts = c(0.8, 0.09, 0.35, 0.6, 0.85),
  scale_factors = c(25, 75, 75, 75, 75)
)

functions_list <- c(fun_1, fun_2, fun_3)

grid <- create_benchmark_grid(
  g_id = c(1, 2, 3),
  n = c(50, 100, 200),
  curvature_percentile = c(15),
  n_lambda = c(10, 20),
  parallel = c(TRUE),
  lambda_grid_spacing = c("log", "sqrt"),
  noise_to_signal = c(0.01, 0.075),
  sigma_amplitude = c(0.2),
  warping = list(
    simple = list(),
    flexible = list(
      warping_gamma = c(3),
      warping_points = c(5)
    )
  ),
  lambda_search_threshold = c(0.1),
  lambda_search_min_bound = c(0.1)
)


benchmark_study <- run_benchmark_study(grid, functions_list,
                                       normalize_functions = TRUE,
                                       output_dir = "../runtime_study",
                                       save_plots = TRUE, runs_per_config = 3,
                                       seed = 1,
                                       parallel = FALSE)
