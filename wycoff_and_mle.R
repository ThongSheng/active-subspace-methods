library(mvtnorm)
library(rstan)
array_id <- as.integer(Sys.getenv('THISJOBVALUE'))
grid <- expand.grid('p' = c(2, 5, 10, 15, 20), 'type' = c('full_rank', '2d'), 
                    'n' = c(25, 75, 125, 175),
                    'seed' = 1:20)

p <- grid$p[array_id] # input dimension
set.seed(grid$seed[array_id])

n <- grid$n[array_id]
x_obs <- matrix(ncol = p, nrow = n, runif(p*n))
plot(x_obs)

dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) -
    matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
dist_tensor_mat_reduced <- dist_tensor_mat[
  upper.tri(matrix(nrow =n, ncol = n), diag = T), ]

compute_cov_mat_from_C <- function(C, dist_tensor_mat_reduced, n,
                                   diagonal_add = .00001) {
  test_distances <- rowSums((dist_tensor_mat_reduced %*% C) *
                              dist_tensor_mat_reduced)
  dist_mat_C <- matrix(nrow = n, ncol = n)
  dist_mat_C[upper.tri(dist_mat_C, diag = T)] <- test_distances
  cov_mat <- exp(-dist_mat_C/2)
  diag(cov_mat) <- diag(cov_mat) + diagonal_add
  cov_mat
}

# We make the eigenvectors randomly; eigenvalues depend on dimension p
W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
if (grid$type[array_id] == 'full_rank') { # all eigenvalvues nonzero
  C <- 500/p * (W_random %*% diag(exp((p:1 -p))) %*% t(W_random))
} else { # only two nonzero eigenvalues
  C <- 500/p * (W_random %*% diag(exp(c(0, -1, rep(-Inf, p-2)))) %*% t(W_random))
}
Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n)
Cov_mat[lower.tri(Cov_mat)] <- t(Cov_mat)[lower.tri(Cov_mat)]
y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))

# start of Wycoff code: uses their activegp package that should be installed
# Wycoff has an input matrix that can sort of like a prior
  # it should probably be set similarly to the prior matrices used for the Bayesian 
prior_scale_matrix <- diag(x = 300/p, nrow = p)

Cov_mat_use <- compute_cov_mat_from_C(prior_scale_matrix, dist_tensor_mat_reduced, n)
Cov_mat_use[lower.tri(Cov_mat_use)] <- t(Cov_mat_use)[lower.tri(Cov_mat_use)]

K_inv <- solve(Cov_mat_use)
K_inv_y <- K_inv %*% y
C_est_wycoff <- matrix(0, p, p)
for (j in 1:p) {
  for (k in 1:p) {
    W_active_gp <- activegp::W_kappa_ij(x_obs, sqrt(1/diag(prior_scale_matrix)), j-1, k-1, 1)
    
    C_est_wycoff[j,k] <- prior_scale_matrix[j,k] - sum(K_inv * W_active_gp) +
      as.vector(t(K_inv_y) %*% W_active_gp %*% K_inv_y)
  }
}
C_est_wycoff

### Start of MLE for general input dimension when observing y

# implement the likelihood function: data | C, which is multivariate normal
likelihood_C <- function(C, dist_tensor_mat_reduced, n,
                         diagonal_add, compute_gradient = F, upper_tri_n = NULL, lower_tri_n = NULL,
                         upper_tri_p = NULL, lower_tri_p = NULL) {
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n,
                                    diagonal_add)
  # The following is the equivalent of mvtnorm::dmvnorm
  # but we do it manually because we only use the upper triangle of Cov_mat
  chol_Cov_mat <- chol(Cov_mat)
  det_Cov_mat <- 2*sum(log(diag(chol_Cov_mat)))
  
  #quadratic form y^T Sigma^-1 y = y^T (L L^T)^-1 y = y^T L^(-T) L^(-1) y 
  backsolved_val <- backsolve(r = chol_Cov_mat, x = y,
                              upper.tri = T, transpose = T)
  cross_Cov_mat <- sum(backsolved_val^2)
  
  # multivariate normal log likelihood:
  -1/2 * ( n * log(2*pi) + det_Cov_mat + cross_Cov_mat)
}

make_C_from_par <- function(par, p) {
  diag_entries_index <- cumsum(1:p)
  diagonal_C <- exp(par[diag_entries_index])
  C_mat_use <- diag(p)
  C_mat_use[upper.tri(C_mat_use)] <- par[-diag_entries_index]
  C_mat_use[lower.tri(C_mat_use)] <- t(C_mat_use)[lower.tri(C_mat_use)]
  diag(sqrt(diagonal_C)) %*% C_mat_use %*% diag(sqrt(diagonal_C)) 
} 

# an example when p = 3
# will get messier for larger p
likelihood_function <- function(par, dist_tensor_mat_reduced, n,
                                diagonal_add) {
  C_mat_use <- make_C_from_par(par, ncol(dist_tensor_mat_reduced))
  #print(C_mat_use)
  if (sum(is.na(C_mat_use)) > 0) {
    return(10^6)
  }
  if (min(eigen(C_mat_use)$value) < 10^-6) {
    return(10^6)
  }
  ll_val <- likelihood_C(C_mat_use, dist_tensor_mat_reduced, n,
                         diagonal_add)
  #print(ll_val)
  -ll_val
}
# initial parameters
init_matrix <- matrix(nrow = p, ncol = p, 1)
diag(init_matrix) <- 300/p
init_par <- log(init_matrix[upper.tri(init_matrix, diag = T)])
make_C_from_par(init_par, p)

# lower and upper bounds
upper_mat <-  matrix(nrow = p, ncol = p, exp(.999))
diag(upper_mat) <- NA
upper_par <- log(upper_mat[upper.tri(upper_mat, diag = T)])

lower_par <- - upper_par

# numerically optimize likelihood, including on the log scale for variance parameters
likelihood_optim <- optim(init_par, likelihood_function,
                          dist_tensor_mat_reduced = dist_tensor_mat_reduced,
                          n = n, method = 'L-BFGS-B',
                          upper = upper_par,
                          lower = lower_par,
                          diagonal_add = .00001)
par <- likelihood_optim$par
likelihood_est <- make_C_from_par(par, p)
likelihood_est
