
library(mvtnorm)
library(rstan)

# test with inverse-wishart multivariate normal

p <- 2 # input dimension
# example C
C <- matrix(nrow = p, ncol = p, c(3, -1.2, -1.2, .5))*20
C <- matrix(nrow = p, ncol = p, c(3, -1, -1, 2.5))*20
C
set.seed(10)

n <- 500

y <-  rmvnorm(n = n, sigma = C)

source('stan/stan_edited_model.R')
m_ss.lniw_mvn   <- stan_model(model_code=sim.sslniw)

it <- 2000 # number of iterations
w <- 500 # number of burn in
C
data_input <- list(y = y, N = n, R = diag(x = 1, nrow = p), # prior
                   k=p, mu0 = rep(0, p),
                   prior_lgn_mean = rep(log(1), p), 
                   prior_lgn_var = rep(log(100), p),
                   prior_dof = p + 3)
out_ln_iw_mvn <- sampling(object=m_ss.lniw_mvn, data = data_input,
                          pars= c('Sigma'), iter = it, chains = 3, warmup=w)
print(out_ln_iw_mvn)
C
plot(out_ln_iw_mvn)
pairs(out_ln_iw_mvn)


# test prior with rescaling by K with inverse-wishart multivariate normal
m_ss.lniw_mvn_rescale   <- stan_model(model_code=sim.sslniw_rescale)
data_input <- list(y = y, N = n, R = diag(x = 1, nrow = p), # prior
                   k=p, mu0 = rep(0, p),
                   prior_lgn_mean = as.array(rep(log(1), p-1)), 
                   prior_lgn_var = as.array(rep(log(100), p-1)),
                   prior_rescale_mean = log(1), 
                   prior_rescale_var = log(2),
                   prior_dof = p + 3)
out_ln_iw_mvn <- sampling(object=m_ss.lniw_mvn_rescale, data = data_input,
                          pars= c('Sigma', 'K'), iter = it, chains = 3, warmup=w)
print(out_ln_iw_mvn)
C
plot(out_ln_iw_mvn)
pairs(out_ln_iw_mvn)

plot(out_ln_iw_mvn, plotfun = 'trace')
plot(out_ln_iw_mvn, plotfun = 'hist')
plot(out_ln_iw_mvn, plotfun = 'rhat')
plot(out_ln_iw_mvn, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(out_ln_iw_mvn@sim$samples[[1]]$`Sigma[1,1]`)



# now try with Gaussian process versions
# sampled points
n <- 500
x_obs <- matrix(ncol = p, nrow = n, runif(p*n))
plot(x_obs)
diag_add = .00001

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
diag(Cov_mat) <- diag(Cov_mat) + diag_add
set.seed(10)
# sample data from multivariate normal
# we are assuming the squared-exponential is the "correct" model
y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))

library(ggplot2)
ggplot(data = data.frame(x_obs, y ), aes(x = X1, X2, color = y)) +
  scale_color_gradient2() + 
  geom_point()

# try without rescaling
m_ss.lniw_gp   <- stan_model(model_code=sim.sslniw_gp)
n_subset <- 100 # only use this many data points
C
data_input_gp <- list(y = y[1:n_subset], N = n_subset, R = diag(x = 1, nrow = 2), # prior
                      k=2, mu0 = rep(0,n_subset), locs = x_obs[1:n_subset,],
                      prior_lgn_mean = rep(log(1), p), 
                      prior_lgn_var = rep(log(100), p),
                      prior_dof = p + 3, diag_add = diag_add)
out_lniw_gp <- sampling(object=m_ss.lniw_gp, data = data_input_gp,
                        pars= c('Sigma'), iter = it, chains = 3, warmup=w)
out_lniw_gp
plot(out_lniw_gp)
pairs(out_lniw_gp)
plot(out_lniw_gp, plotfun = 'trace')
plot(out_lniw_gp, plotfun = 'hist')
plot(out_lniw_gp, plotfun = 'rhat')
plot(out_lniw_gp, show_density = TRUE, ci_level = 0.5, fill_color = "purple")
plot(out_lniw_gp@sim$samples[[1]]$`Sigma[1,1]`)


# try with rescaling
m_ss.lniw_gp_rescale   <- stan_model(model_code=sim.sslniw_gp_rescale)
n_subset <- 100
ggplot(data = data.frame(x_obs, y )[1:n_subset,], aes(x = X1, X2, color = y)) +
  scale_color_gradient2() + 
  geom_point()
data_input_gp_rescale <- list(y = y[1:n_subset], N = n_subset, R = diag(x = 1, nrow = 2), # prior
                              k=2, mu0 = rep(0,n_subset), locs = x_obs[1:n_subset,],
                              prior_lgn_mean = array(rep(log(.5), p-1)), 
                              prior_lgn_var = array(rep(log(100), p-1)),
                              prior_dof = p + 3,
                              prior_rescale_mean = log(1), 
                              prior_rescale_var = log(2), diag_add = .00001)
cor_inits<- runif(3, -.5, .5) # the third chain often starts with really high correlation, leading to bad mixing
out_lniw_gp_rescale <- sampling(object=m_ss.lniw_gp_rescale, data = data_input_gp_rescale,
                                pars= c('Sigma'), iter = it, chains = 3, warmup=w, 
                                init = list(list('xi' = array(.5), 
                                                 'Q1' = matrix(c(1,cor_inits[1], cor_inits[1], 1), nrow = 2)),
                                            list('xi' = array(.5), 
                                                 'Q1' = matrix(c(1,cor_inits[2], cor_inits[2], 1), nrow = 2)),
                                            list('xi' = array(.5), 
                                                 'Q1' = matrix(c(1,cor_inits[3], cor_inits[3], 1), nrow = 2))))
print(out_lniw_gp_rescale)
plot(out_lniw_gp_rescale)
pairs(out_lniw_gp_rescale)
out_lniw_gp_rescale@inits[[1]]
out_lniw_gp_rescale@inits[[2]]
out_lniw_gp_rescale@inits[[3]]

# compare chains
plot(out_lniw_gp_rescale@sim$samples[[1]]$`Sigma[1,1]`)
abline(h = C[1,1], col = 2)
plot(out_lniw_gp_rescale@sim$samples[[2]]$`Sigma[1,1]`)
abline(h = C[1,1], col = 2)
plot(out_lniw_gp_rescale@sim$samples[[3]]$`Sigma[1,1]`)
abline(h = C[1,1], col = 2)
plot(out_lniw_gp_rescale@sim$samples[[1]]$`Sigma[1,2]`)
abline(h = C[1,2], col = 2)
plot(out_lniw_gp_rescale@sim$samples[[2]]$`Sigma[1,2]`)
abline(h = C[1,2], col = 2)
plot(out_lniw_gp_rescale@sim$samples[[3]]$`Sigma[1,2]`)
abline(h = C[1,2], col = 2)
plot(out_lniw_gp_rescale@sim$samples[[1]]$`Sigma[2,2]`)
abline(h = C[2,2], col = 2)
plot(out_lniw_gp_rescale@sim$samples[[2]]$`Sigma[2,2]`)
abline(h = C[2,2], col = 2)
plot(out_lniw_gp_rescale@sim$samples[[3]]$`Sigma[2,2]`)
abline(h = C[2,2], col = 2)

