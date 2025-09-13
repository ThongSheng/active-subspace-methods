

# Gautier
library(mvtnorm)

# --- Simulation Grid Definition ---
grid <- expand.grid(
  'p' = c(2, 10, 20),
  'n' = c(20, 100, 150),
  'seed' = 1:3,
  'func' = c('determ_full', 'determ_2d', 'GP_full', 'GP_2d'))

# plot(seq(0, 20, length.out = 400), dlnorm(seq(0, 20, length.out = 400), 
#                                           meanlog = 0, sdlog = sqrt(log(2))), type = 'l')
# plot(seq(0, 20, length.out = 400), dlnorm(seq(0, 20, length.out = 400), 
#                                           meanlog = 0, sdlog = sqrt(log(2))), type = 'l')# Extract parameters for this specific job array run
array_id <- as.integer(Sys.getenv('THISJOBVALUE'))
array_id <- 1
p <- grid$p[array_id]
n <- grid$n[array_id]
n <- 20
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

if (grid$func[array_id] == 'determ_full') {
  f <- function(x) {sum(x^2)/sqrt(p)}
  C <- matrix(0, nrow = p, ncol = p)
  diag(C) <- 4/(3*p)
  if (p > 1) {
    C[upper.tri(C) | lower.tri(C)] <- 1/p
  }
  y <- apply(x_obs, 1, f)
} else if (grid$func[array_id] == 'determ_2d') {
  f <- function(x) {sum(x[1:2]^2)/sqrt(p)}
  C <- matrix(0, nrow = p, ncol = p)
  if (p >= 2) {
    C[1,1] <- C[2,2] <- 4/(3*p)
    C[1,2] <- C[2,1] <- 1/p
  }
  y <- apply(x_obs, 1, f)  
} else if (grid$func[array_id] == 'GP_full') {
  W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
  C <- 500/p * (W_random %*% diag(exp((p:1 -p))) %*% t(W_random))
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n)
  Cov_mat[lower.tri(Cov_mat)] <- t(Cov_mat)[lower.tri(Cov_mat)]
  y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))
  
} else if (grid$func[array_id] == 'GP_2d') {
  W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
  C <- 500/p * (W_random %*% diag(exp(c(0, -1, rep(-Inf, p-2)))) %*% t(W_random))
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n)
  Cov_mat[lower.tri(Cov_mat)] <- t(Cov_mat)[lower.tri(Cov_mat)]
  y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))  
}
library(rstan)
stan_model_object <- stan_model(model_code = sim.gp_reduce)
#save(stan_model_object, file = 'test_stan_model.RData')
# MCMC initializations
n_chains <- 4
it <- 700
w <- 300

# Run STAN sampling
a_time <- Sys.time()
out_model <- sampling(object = stan_model_object, data = list('N' = n, 'k' = p, 'k_reduce' = 2, 'y' = y,
                             'mu0' = rep(0, n), 'locs' = x_obs, 'diag_add' = .00001),
                      pars = c('Sigma'), iter = it, chains = n_chains, warmup = w,
                      #init = initial_values,
                      cores = n_chains)

#pairs(out_model)
eigen(C)$vectors[,1]

C_est_sample <- extract(out_model)$Sigma

post_mean <- apply(C_est_sample, c(2,3), mean)
post_mean
C

ev1_est <- eigen(post_mean)$vectors[,1]
ev2_est <- eigen(post_mean)$vectors[,2]

ev1 <- eigen(C)$vectors[,1]
ev2 <- eigen(C)$vectors[,2]
similarity1 <- abs(sum(ev1_est * ev1))
similarity2 <- abs(sum(ev2_est * ev2))

