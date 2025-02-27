source("./R/sim_data.R")

t_grid <- seq(0, 1, length.out = 1000)
fun_1 <- generate_benchmark_function(list("normal", "normal"),
                                     c(10, 10), c(0.3, 0.7), c(10, 10))

fun_1_tf <- tf::tfd(fun_1(t_grid), t_grid)
plot(fun_1_tf, t_grid)

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

fun_2_tf <- tf::tfd(fun_2(t_grid), t_grid)
plot(fun_2_tf, t_grid)


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

fun_3_tf <- tf::tfd(fun_3(t_grid), t_grid)
plot(fun_3_tf, t_grid)
