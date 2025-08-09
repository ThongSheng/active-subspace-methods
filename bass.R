library(mvtnorm)
library(BASS)
library(concordance)

# --- Simulation Grid Definition ---
grid <- expand.grid(
  'p' = c(2, 10, 20),
  'type' = c('full_rank', '2d'),
  'n' = c(20, 100, 150),
  'seed' = 1:3
)

# Extract parameters for this specific job array run
array_id <- as.integer(Sys.getenv('THISJOBVALUE'))
p <- grid$p[array_id]
n <- grid$n[array_id]
set.seed(grid$seed[array_id])

# --- Data Generation  ---
x_obs <- matrix(ncol = p, nrow = n, runif(p*n))
dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) - matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
dist_tensor_mat_reduced <- dist_tensor_mat[upper.tri(matrix(nrow =n, ncol = n), diag = T), ]

compute_cov_mat_from_C <- function(C, dist_tensor_mat_reduced, n, diagonal_add = .00001) {
  test_distances <- rowSums((dist_tensor_mat_reduced %*% C) * dist_tensor_mat_reduced)
  dist_mat_C <- matrix(nrow = n, ncol = n)
  dist_mat_C[upper.tri(dist_mat_C, diag = T)] <- test_distances
  cov_mat <- exp(-dist_mat_C/2)
  diag(cov_mat) <- diag(cov_mat) + diagonal_add
  cov_mat
}

W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
if (grid$type[array_id] == 'full_rank') {
  C <- 500/p * (W_random %*% diag(exp((p:1 -p))) %*% t(W_random))
} else {
  C <- 500/p * (W_random %*% diag(exp(c(0, -1, rep(-Inf, p-2)))) %*% t(W_random))
}
Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n)
Cov_mat[lower.tri(Cov_mat)] <- t(Cov_mat)[lower.tri(Cov_mat)]
y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))

# --- Sampling ---
a_time <- Sys.time()
mod_bass <- BASS::bass(x_obs, y, verbose=FALSE)
pred_C <- C_bass(mod_bass)
b_time <- Sys.time()
time_used <- b_time - a_time
units(time_used) <- 'mins'

# Save results
output_dir <- '/scratch/negishi/angt/simulation_results/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
save_file_name <- paste0(output_dir, "BASS_", array_id, '.RData')
save(pred_C, C, y, time_used, file = save_file_name)
q(save = 'no')
