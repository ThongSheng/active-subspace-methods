
p <- 3 # input dimension
# example C
C <- matrix(nrow = p, ncol = p, c(1, -.4, .2, -.4, 3, 1.8,.2, 1.8, 6))*25
C
set.seed(10)
# sampled points
n <- 200
x_obs <- matrix(ncol = p, nrow = n, runif(p*n))
plot(x_obs)

# We will reorganize the distances theta - theta'
dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) -
    matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}

t(dist_tensor[1,2,]) %*% C %*% dist_tensor[1,2,]

# reshape from n times n times p; to n^2 times p
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)

# remove lower diagonal to not repeat calculations
dist_tensor_mat_reduced <- dist_tensor_mat[
  upper.tri(matrix(nrow =n, ncol = n), diag = T), ]

compute_cov_mat_from_C <- function(C, dist_tensor_mat_reduced, n,
                                   diagonal_add = .00001) {
  # what we want is the diagonal of 
  # dist_tensor_mat_reduced %*% C %*% t(dist_tensor_mat_reduced)
  # the rowsums approach computes only the diagonal of this matrix
  test_distances <- rowSums((dist_tensor_mat_reduced %*% C) *
                              dist_tensor_mat_reduced)
  # we put this back into n times n matrix form
  dist_mat_C <- matrix(nrow = n, ncol = n)
  # we only compute the upper diagonal
  dist_mat_C[upper.tri(dist_mat_C, diag = T)] <- test_distances
  # add a bit to diagonal for stability
  cov_mat <- exp(- dist_mat_C/2)
  diag(cov_mat) <- diag(cov_mat) + diagonal_add
  cov_mat
}
Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n)


library(mvtnorm)
set.seed(10)
# sample data from multivariate normal
# we are assuming the squared-exponential is the "correct" model
# use cholesky to only use upper triangle
chol_Cov_mat <- chol(Cov_mat)
#y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))
y <- t(chol_Cov_mat) %*% rnorm(n)

library(ggplot2)

ggplot(data = data.frame(x_obs, y)) +
  geom_point(aes(x = X1, X2, color = y)) +
  scale_color_gradient2()
ggplot(data = data.frame(x_obs, y)) +
  geom_point(aes(x = X1, X3, color = y)) +
  scale_color_gradient2()
ggplot(data = data.frame(x_obs, y)) +
  geom_point(aes(x = X3, X2, color = y)) +
  scale_color_gradient2()


# implement the likelihood function: data | C, which is multivariate normal
likelihood_C <- function(C, dist_tensor_mat_reduced, n,
                         diagonal_add) {
  Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n,
                                    diagonal_add)
  # The following is the equivalent of mvtnorm::dmvnorm
  # but we do it manually because we only use the upper triangle of Cov_mat
  chol_Cov_mat <- chol(Cov_mat)
  
  # log determinant, log(det(Sigma)) = log( det(Sigmachol)^2 )
  #                    = 2 * log( det(Sigmachol))
  #                    = 2 * log(prod(diag(Sigmachol))) # triangular
  #                    = 2 * sum(log(diag(Sigmachol))))
  det_Cov_mat <- 2*sum(log(diag(chol_Cov_mat)))
  
  #quadratic form y^T Sigma^-1 y = y^T (L L^T)^-1 y = y^T L^(-T) L^(-1) y 
  cross_Cov_mat <- sum(backsolve(r = chol_Cov_mat, x = y,
                                 upper.tri = T, transpose = T)^2)
  
  # multivariate normal log likelihood:
  -1/2 * ( n * log(2*pi) + det_Cov_mat + cross_Cov_mat)
}

# an example when p = 3
# will get messier for larger p
likelihood_function <- function(par, dist_tensor_mat_reduced, n,
                                diagonal_add) {
  diagonal_C <- exp(par[c(1,3,6)])
  C_mat_use <- matrix(nrow = 3, ncol =3, 
                      c(diagonal_C[1], 
                        sqrt(diagonal_C[1] * diagonal_C[2]) *par[2],
                        sqrt(diagonal_C[1] * diagonal_C[3]) *par[4],
                        sqrt(diagonal_C[1] * diagonal_C[2]) *par[2],
                        diagonal_C[2], 
                        sqrt(diagonal_C[2] * diagonal_C[3]) *par[5],
                        sqrt(diagonal_C[1] * diagonal_C[3]) *par[4],
                        sqrt(diagonal_C[2] * diagonal_C[3]) *par[5],
                        diagonal_C[3]))
  if (min(eigen(C_mat_use)$value) < 10^-6) {
    return(10^6)
  }
  print(C_mat_use)
  ll_val <- likelihood_C(C_mat_use, dist_tensor_mat_reduced, n,
                         diagonal_add)
  print(ll_val)
  -ll_val
}
# numerically optimize likelihood, including on the log scale for variance parameters
likelihood_optim <- optim(c(log(120), (.2), log(50), .2, .2, log(50)), likelihood_function,
                          dist_tensor_mat_reduced = dist_tensor_mat_reduced,
                          n = n, method = 'L-BFGS-B',
                          upper = c(NA, .999, NA, .999, .999, NA),
                          lower = c(NA, -.999, NA, -.999, -.999, NA),
                          diagonal_add = .00001)
par <- likelihood_optim$par
diagonal_C <- exp(par[c(1,3,6)])
likelihood_est <- matrix(nrow = 3, ncol =3, 
                         c(diagonal_C[1], 
                           sqrt(diagonal_C[1] * diagonal_C[2]) *par[2],
                           sqrt(diagonal_C[1] * diagonal_C[3]) *par[4],
                           sqrt(diagonal_C[1] * diagonal_C[2]) *par[2],
                           diagonal_C[2], 
                           sqrt(diagonal_C[2] * diagonal_C[3]) *par[5],
                           sqrt(diagonal_C[1] * diagonal_C[3]) *par[4],
                           sqrt(diagonal_C[2] * diagonal_C[3]) *par[5],
                           diagonal_C[3]))


likelihood_est
C
likelihood_est/C
# we have direct likelihood estimation in the squared exponential covariance and it works well

# Bayesian MCMC using this
# metropolis hasting proposes a new C matrix, then accepts/rejects this proposal
# we use a Wishart distribution for the proposal



# set prior distribution
prior_C <- matrix(nrow = p, ncol = p , c(30, 20, 20, 
                                         20, 80, 20, 
                                         20, 20, 160))
#prior_C <- C
df_prior <- p + 10 # degrees of freedom

MC <- 
  5000 # number of Monte carlo iterations
C_MCMC <- array(dim = c(p, p, MC))
C_MCMC[,,1] <- diag(100, nrow = p)
C_MCMC[,,1] <- diag(100 + rnorm(p), nrow = p)
#C_MCMC[,,1] <- C # initial start
accept_prob <- rep(0, MC)
accept <- rep(0, MC)

# prior and likelihood values at initialization
prior <- 'ShrinkageInvWishart'
# prior and likelihood values at initialization
if (prior == 'InvWishart') {
  prior_prev <- CholWishart::dInvWishart(C_MCMC[,,1] / df_prior, 
                                         Sigma = prior_C, df = df_prior,
                                         log = T)
} else if (prior == 'ShrinkageInvWishart') {
  eigen_C <- eigen(C_MCMC[,,1] / df_prior)$values
  eigen_C_mat <- matrix(eigen_C, nrow = p, ncol = p) - matrix(eigen_C, nrow = p, ncol = p, byrow = T)
  upper_tri_part <- upper.tri(eigen_C_mat, diag = F)
  prior_prev <- CholWishart::dInvWishart(C_MCMC[,,1] / df_prior, 
                                         Sigma = prior_C, df = df_prior,
                                         log = T) + 
    -sum(log(eigen_C_mat[upper_tri_part]))
}
likelihood_value_prev <- likelihood_C(C_MCMC[,,1], 
                                      dist_tensor_mat_reduced, n,
                                      diagonal_add = .00001)

df_prop <- p + 1000 # degrees of freedom for proposal distribution

for (mc in 2:MC) {
  
  # propose next sample
  C_prop <- rWishart(n = 1, df = df_prop, C_MCMC[,,mc - 1])[,,1]/df_prop
  
  # compute log-likelihood/prior for proposed C
  likelihood_value_prop <- likelihood_C(C_prop, dist_tensor_mat_reduced, n,
                                        diagonal_add = .00001)
  if (prior == 'InvWishart') {
    prior_prop <- CholWishart::dInvWishart(C_prop / df_prior, 
                                           Sigma = prior_C, df = df_prior,
                                           log = T)
  } else if (prior == 'ShrinkageInvWishart') {
    eigen_C <- eigen(C_prop / df_prior)$values
    prior_prop <- CholWishart::dInvWishart(C_prop / df_prior, 
                                           Sigma = prior_C, df = df_prior,
                                           log = T) -
      sum(log(eigen_C_mat[upper_tri_part]))
  }
  # adjustment for MCMC
  q_prop_prev <- CholWishart::dWishart(C_prop * df_prop, Sigma = C_MCMC[,,mc - 1], df = df_prop,
                                       log = T)
  q_prev_prop <- CholWishart::dWishart(C_MCMC[,,mc - 1] * df_prop, Sigma = C_prop, df = df_prop,
                                       log = T)
  
  accept_prob[mc] <- exp(likelihood_value_prop - likelihood_value_prev + 
                           prior_prop - prior_prev +
                           q_prev_prop - q_prop_prev)
  accept_prob[mc] <- ifelse(accept_prob[mc] > 1, 1, accept_prob[mc])
  if (runif(1) < accept_prob[mc]) {
    C_MCMC[,,mc] <- C_prop
    accept[mc] <- 1
    likelihood_value_prev <- likelihood_value_prop
    prior_prev <- prior_prop
  } else {
    C_MCMC[,,mc] <- C_MCMC[,,mc - 1] 
  }
  print(likelihood_value_prop + prior_prop)
  print(C_MCMC[,,mc])
}

plot(C_MCMC[1,1,], cex = .2)
abline(h = C[1,1], col = 2)
plot(C_MCMC[2,1,], cex = .2)
abline(h = C[2,1], col = 2)
plot(C_MCMC[3,1,], cex = .2)
abline(h = C[3,1], col = 2)
plot(C_MCMC[2,2,], cex = .2)
abline(h = C[2,2], col = 2)
plot(C_MCMC[2,3,], cex = .2)
abline(h = C[2,3], col = 2)
plot(C_MCMC[3,3,], cex = .2)
abline(h = C[3,3], col = 2)

# remove first samples due to initialization
burn_in_number <- 1000
C_MCMC_with_burn_in <- C_MCMC
C_MCMC <- C_MCMC[,,(burn_in_number+1):(dim(C_MCMC)[3])]

hist(accept_prob)
mean(accept)
apply(C_MCMC, c(1,2), mean)
apply(C_MCMC, c(1,2), var)
plot(C_MCMC[1,1,], cex = .2)
abline(h = C[1,1], col = 2)
plot(C_MCMC[2,1,], cex = .2)
abline(h = C[2,1], col = 2)
plot(C_MCMC[2,2,], cex = .2)
abline(h = C[2,2], col = 2)

hist(C_MCMC[1,1,], cex = .2)
abline(v = C[1,1], col = 2)
hist(C_MCMC[1,2,], cex = .2)
abline(v = C[1,2], col = 2)
hist(C_MCMC[2,2,], cex = .2)
abline(v = C[2,2], col = 2)

# check eigenvalues
eig_values <- apply(C_MCMC, 3, function(x) eigen(x)$values)
hist(eig_values[1,])
abline(v = eigen(C)$values[1], col = 2)
hist(eig_values[2,])
abline(v = eigen(C)$values[2], col = 2)
hist(eig_values[3,])
abline(v = eigen(C)$values[3], col = 2)

# cosine similarity of eigenvectors
vectors_cos_sim <- apply(C_MCMC, 3, function(x) sum(eigen(x)$vectors[,1] * eigen(C)$vectors[,1])/
                           sqrt(sum(eigen(x)$vectors[,1]^2) * sum(eigen(C)$vectors[,1]^2)))
hist(vectors_cos_sim)
hist(abs(vectors_cos_sim))

# second eigenvector
vectors_cos_sim <- apply(C_MCMC, 3, function(x) sum(eigen(x)$vectors[,2] * eigen(C)$vectors[,2])/
                           sqrt(sum(eigen(x)$vectors[,2]^2) * sum(eigen(C)$vectors[,2]^2)))
hist(vectors_cos_sim)
hist(abs(vectors_cos_sim))


