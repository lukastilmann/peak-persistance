library(R.matlab)
library(tf)
library(tidyfun)
devtools::load_all(".")

mat_data <- readMat("./data/example1.mat")
t_grid <- as.vector(mat_data$t)
ground_truth <- mat_data$g
ground_truth_tf <- tf::tfd(ground_truth, t_grid)
unaligned_curves <- mat_data$f
unaligned_curves_tf <- tf::tfd(unaligned_curves, t_grid)

plot(unaligned_curves_tf, t_grid)
ppd_result <- peak_persistance_diagram(unaligned_curves_tf, t_grid, lambda_search_threshold = 0.15)

mfn <- ppd_result$mean_function
peak_locs <- ppd_result$peak_locs
valley_locs <- ppd_result$valley_locs
significant_peaks <- ppd_result$significant_peaks
peak_labels <- ppd_result$labels
t_grid <- ppd_result$time_grid
fun_aligned <- ppd_result$aligned_functions

k = 8
rho = 1e-11
basis_type = "tp"

fn_est <- shape_constrained_estimation(fun_aligned, peak_locs, valley_locs,
                                       significant_peaks, peak_labels, mfn, k,
                                       t_grid, rho, basis_type)

vis_df <- data.frame(fun_aligned)

ggplot(vis_df, aes(y = fun_aligned)) +
  geom_spaghetti(alpha = 0.3) +
  geom_spaghetti(aes(y = fn_est), color = "orange", linewidth = 1.5) +
  geom_spaghetti(aes(y = mfn), color = "green", linewidth = 1.5) +
  geom_spaghetti(aes(y = ground_truth_tf), color = "red", linewidth = 1.5)
