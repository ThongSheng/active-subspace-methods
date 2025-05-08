
p <- 2 # input dimension
# example C
C <- matrix(nrow = p, ncol = p, c(3, .2, .2, .5))*50
C
set.seed(10)
# sampled points
gridded <- F 
if (gridded) {
  grid_vals <- seq(0, 1, length.out = 35)
  sqrt_n <- length(grid_vals)
  n = sqrt_n^2
  x_obs <- as.matrix(expand.grid(grid_vals, grid_vals))
} else {
  n <- 200
  x_obs <- matrix(ncol = p, nrow = n, runif(p*n))
}
plot(x_obs)

# make matrices of distances in the direction of each dimension
unscaled_dist1 <- matrix(x_obs[,1], nrow = n, ncol = n) - matrix(x_obs[,1], nrow = n, ncol = n, byrow = T)
unscaled_dist2 <- matrix(x_obs[,2], nrow = n, ncol = n) - matrix(x_obs[,2], nrow = n, ncol = n, byrow = T)

# create covariance matrix with squared exponential kernel
  # and anisotropy matrix C
dist_sq_mat <- unscaled_dist1^2 *C[1,1] + 
  unscaled_dist2^2 *C[2,2] + 
  unscaled_dist2 * unscaled_dist1 *C[1,2] +
  unscaled_dist2 * unscaled_dist1 *C[2,1]

Cov_mat <- exp(- dist_sq_mat/2)


# add a bit to diagonal for stability
diag(Cov_mat) <- diag(Cov_mat) + .00001

library(mvtnorm)
set.seed(10)
# sample data from multivariate normal
  # we are assuming the squared-exponential is the "correct" model
y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))


# implement the likelihood function: data | C, which is multivariate normal
likelihood_C <- function(C) {
  # for d > 2, we will need to generalize; will need efficient matrix implementation
  dist_sq_mat <- unscaled_dist1^2 * C[1,1] + 
    unscaled_dist2^2 * C[2,2] + 
    unscaled_dist2 * unscaled_dist1 * C[1,2] +
    unscaled_dist2 * unscaled_dist1 * C[2,1]
  Cov_mat <- exp(-dist_sq_mat/2)
  diag(Cov_mat) <- diag(Cov_mat) + .00001
  mvtnorm::dmvnorm(y, sigma = Cov_mat, log = T)
}

likelihood_function <- function(par) {
  C_mat_use <- matrix(nrow = 2, ncol =2, c(exp(par[1]), par[2], par[2], exp(par[3])))
  print(C_mat_use)
  ll_val <- likelihood_C(C_mat_use)
  print(ll_val)
  -ll_val
}
# numerically optimize likelihood, including on the log scale for variance parameters
likelihood_optim <- optim(c(log(120), (20), log(50)), likelihood_function)

likelihood_est <- matrix(nrow = p, ncol =p, c(exp(likelihood_optim$par[1]), likelihood_optim$par[2],
                                              likelihood_optim$par[2], exp(likelihood_optim$par[3])))
likelihood_est
C
likelihood_est/C
mvtnorm::dmvnorm(y, sigma = Cov_mat, log = T)

# we have direct likelihood estimation in the squared exponential covariance and it works well

# Bayesian MCMC using this
  # metropolis hasting proposes a new C matrix, then accepts/rejects this proposal
    # we use a Wishart distribution for the proposal



# set prior distribution
prior_C <- matrix(nrow = p, ncol = p , c(200, 20, 20, 20))
#prior_C <- C
df_prior <- p + 10 # degrees of freedom

MC <- 
  5000 # number of Monte carlo iterations
C_MCMC <- array(dim = c(p, p, MC))
C_MCMC[,,1] <- c(1, 0, 0, 1)
#C_MCMC[,,1] <- C # initial start
accept_prob <- rep(0, MC)
accept <- rep(0, MC)

# prior and likelihood values at initialization
prior_prev <- CholWishart::dInvWishart(C_MCMC[,,1] / df_prior, 
                                       Sigma = prior_C, df = df_prior,
                                       log = T)
likelihood_value_prev <- likelihood_C(C_MCMC[,,1])

df_prop <- p + 1000 # degrees of freedom for proposal distribution

for (mc in 2:MC) {
  
  # propose next sample
  C_prop <- rWishart(n = 1, df = df_prop, C_MCMC[,,mc - 1])[,,1]/df_prop
  
  # compute log-likelihood/prior for proposed C
  likelihood_value_prop <- likelihood_C(C_prop)
  prior_prop <- CholWishart::dInvWishart(C_prop / df_prior, 
                                         Sigma = prior_C, df = df_prior,
                                         log = T)
  
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
plot(C_MCMC[2,2,], cex = .2)
abline(h = C[2,2], col = 2)

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

# cosine similarity of eigenvectors
vectors_cos_sim <- apply(C_MCMC, 3, function(x) sum(eigen(x)$vectors[,1] * eigen(C)$vectors[,1])/
                           sqrt(sum(eigen(x)$vectors[,1]^2) * sum(eigen(C)$vectors[,1]^2)))
hist(vectors_cos_sim)


