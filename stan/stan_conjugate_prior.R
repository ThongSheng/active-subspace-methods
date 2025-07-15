library(mvtnorm)
library(rstan)
array_id <- as.integer(Sys.getenv('THISJOBVALUE'))
grid <- expand.grid('p' = c(2, 10, 15, 20),  
                    'type' = c('full_rank', '2d'),
                    'add_dof' = c(2, 10, 50),
                    'n' = c(20, 120, 220, 320),
                    'seed' = 1:20)

p <- grid$p[array_id] # input dimension
set.seed(grid$seed[array_id])

n <- grid$n[array_id]
x_obs <- matrix(ncol = p, nrow = n, runif(p*n))
plot(x_obs)

dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) - matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
# dist_tensor_mat_reduced <- dist_tensor_mat[
#   upper.tri(matrix(nrow =n, ncol = n), diag = T), ]
dist_tensor_mat_reduced_crossprod <- lapply(1:n, 
                                            function(i) lapply(1:n, function(j) {
                                              tcrossprod(dist_tensor[i,j,])
                                            }))

# compute_cov_mat_from_C <- function(C, dist_tensor_mat_reduced, n, diagonal_add = .00001) {
#   test_distances <- rowSums((dist_tensor_mat_reduced %*% C) *
#                               dist_tensor_mat_reduced)
#   dist_mat_C <- matrix(nrow = n, ncol = n)
#   dist_mat_C[upper.tri(dist_mat_C, diag = T)] <- test_distances
#   cov_mat <- exp(-dist_mat_C/2)
#   diag(cov_mat) <- diag(cov_mat) + diagonal_add
#   cov_mat
# }

# replaced compute_cov_mat_from_C with compute_grad_cov_mat_from_C
compute_grad_cov_mat_from_C <- function(C, dist_tensor_mat, n, diagonal_add = .00001) {
  test_distances <- matrix(rowSums((dist_tensor_mat %*% C) * dist_tensor_mat), n, n)
  exp_mat <- exp(-test_distances/2)
  NA_mat <- matrix(NA, p, p)
  all_mats <- lapply(1:n, 
                     function(i) {lapply(1:n, function(j) {
                       if (j < i) {
                         return(NA_mat) # for the lower triangle, do not compute, saves computation time
                       } else {
                         mat <- - ((C %*% dist_tensor_mat_reduced_crossprod[[i]][[j]] %*% C - C)*exp_mat[i,j])
                         if (i == j) {
                           diag(mat) <- diag(mat) + diagonal_add
                         }
                         return(mat)
                       }})
                     })
  # combine all matrices into one pn times pn matrix
  cov_mat <- do.call(rbind, lapply(all_mats, function(x) do.call(cbind, x)))
}

# We make the eigenvectors randomly; eigenvalues depend on dimension p
W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
if (grid$type[array_id] == 'full_rank') { # all eigenvalvues nonzero
  C <- 500/p * (W_random %*% diag(exp((p:1 -p))) %*% t(W_random))
} else { # only two nonzero eigenvalues
  C <- 500/p * (W_random %*% diag(exp(c(0, -1, rep(-Inf, p-2)))) %*% t(W_random))
}

Cov_mat <- compute_grad_cov_mat_from_C(C, dist_tensor_mat, n)
chol_Cov_mat <- chol(Cov_mat)
y <- t(chol_Cov_mat) %*% rnorm(n * p)
grad <- matrix(y, nrow = n, ncol = p, byrow = TRUE)

source('stan/stan_conjugate_prior_model.R')
if (F) {
  m_ss.conj_prior <- stan_model(model_code=sim.conj_prior)
  save(m_ss.conj_prior, file = 'stan/conjugate_prior_model.RData')
} else {
  load('stan/conjugate_prior_model.RData')
}

n_chains <- 4
it <- 1500 # number of iterations
w <- 500 # number of burn in

data_input <- list(y = grad, 
                   N = n, 
                   prior_scale = diag(x = 1, nrow = p), 
                   k=p, 
                   mu0 = rep(0, p),
                   prior_dof = p + grid$add_dof[array_id])

cor_inits <- runif(n_chains, 0, .5) # the third chain often starts with really high correlation, leading to bad mixing
cor_inits_mat <- lapply(cor_inits, function(x) {
  mat <- matrix(x, nrow = p, ncol = p)
  diag(mat) <- 1
  list('Sigma' = mat)
})
sapply(cor_inits_mat, function(x) min(eigen(x[[1]])$values))
a_time <- Sys.time()
out_conj_prior <- sampling(object=m_ss.conj_prior, data = data_input,
                          pars= c('Sigma'), iter = it, chains = n_chains, warmup=w,
                          init = cor_inits_mat, 
                          cores = n_chains)
extract_vals <- extract(out_conj_prior) # extract samples
summary_vals <- summary(out_conj_prior) # extract posterior means/etc.
b_time <- Sys.time()
time_used <- b_time - a_time
units(time_used) <- 'mins'

save(extract_vals, summary_vals, C, y, time_used,
     file = paste0('/scratch/negishi/angt/init_conj_prior/', array_id, '.RData'))
q(save = 'no')

### after running, process results here

files <- list.files('/scratch/negishi/angt/init_conj_prior')
C_mat <- list()
p_vals <- rep(NA, length(files))
row_number <- as.integer(sapply(strsplit(files, '\\.'), function(x) x[[1]]))
cos_sim <- matrix(nrow = length(files), ncol = 4000) 
cos_sim2 <- matrix(nrow = length(files), ncol = 4000)
cos_sim3 <- matrix(nrow = length(files), ncol = 4000)
cos_sim_ref <- matrix(nrow = length(files), ncol = 4000)
mean_eigen <- matrix(nrow = length(files), ncol = 5)
time_used_vals <- rep(NA, length(files))
posterior_means <- list()
for (i in 1:length(files)) {
  load(paste0('/scratch/negishi/angt/init_conj_prior/', files[i]))
  C_mat[[i]] <- C
  p_vals[i] <- nrow(C)
  #  posterior mean of Sigma
  posterior_means[[i]] <- apply(extract_vals$Sigma, c(2,3), mean) 
  # eigenstructure of sampled C/Sigma matrices
  Sigma_mats <- extract_vals$Sigma
  true_eigen <- eigen(C, symmetric = T)$vectors[,1]
  true_eigen2 <- eigen(C, symmetric = T)$vectors[,2]
  
  eigen_Sigma <- lapply(1:(dim(Sigma_mats)[1]), function(x) eigen(Sigma_mats[x,,], symmetric = T))
  eigen_values <- sapply(eigen_Sigma, function(x) x$values)
  eigen_vectors <- lapply(eigen_Sigma, function(x) x$vectors)
  
  # cosine similarity between first true eigenvector and sampled eigenvectors
  cos_sim[i,] <- sapply(1:length(eigen_vectors), function(x) abs(sum(eigen_vectors[[x]][,1] *
                                                                       true_eigen)))
  cos_sim2[i,] <- sapply(1:length(eigen_vectors), function(x) abs(sum(eigen_vectors[[x]][,2] *
                                                                        true_eigen2)))
  # if p>2, also look at the 3rd eigenvector
  if (p_vals[i] > 2) {
    true_eigen3 <- eigen(C, symmetric = T)$vectors[,3]
    cos_sim3[i,] <- sapply(1:length(eigen_vectors), function(x) abs(sum(eigen_vectors[[x]][,3] *
                                                                          true_eigen3)))
    mean_eigen[i,] <- rowMeans(eigen_values)[1:5]
  } else {
    mean_eigen[i,] <- c(rowMeans(eigen_values), NA, NA, NA)

  }
  # compare with randomly generated directions (a baseline)
    # on average vectors in 20 dimensions are less similar than vectors in 2 dimensions
  cos_sim_ref[i,] <- sapply(1:length(eigen_vectors), 
                            function(x) {test <- rnorm(p_vals[i]);abs(sum(test * true_eigen))/ 
                              sqrt(sum(test^2))})
  time_used_vals[i] <- as.double(time_used)
  print(i)
}


grid <- expand.grid('p' = c(2, 10, 15, 20),  
                    'type' = c('full_rank', '2d'),
                    'add_dof' = c(2, 10, 50),
                    'n' = c(20, 120, 220, 320),
                    'seed' = 1:20)
grid$index <- 1:nrow(grid)

grid_reduced <- grid[row_number,]
grid_reduced$index2 <- 1:nrow(grid_reduced)
960/6
table(grid_reduced$p)
table(grid_reduced$n)
table(grid_reduced$type)
table(grid_reduced$add_dof)
library(ggplot2)

ggplot(data = data.frame(grid_reduced, value = time_used_vals)) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p, labeller = label_both) +
  labs(x = 'Sample size (n)', y = 'Time to compute (minutes)')


ggplot(data = data.frame(grid_reduced, value = rowMeans(cos_sim))) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p, labeller = label_both) +
  geom_hline(data = data.frame(grid_reduced, value = rowMeans(cos_sim_ref)), aes(yintercept = value),
             linewidth = .01)+
  labs(x = 'Sample size (n)', 
       y = 'Cosine similarity of\nfirst active\nsubspace direction',
       color = 'C matrix rank')

ggplot(data = data.frame(grid_reduced, value = rowMeans(cos_sim2))) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p, labeller = label_both) +
  geom_hline(data = data.frame(grid_reduced, value = rowMeans(cos_sim_ref)), aes(yintercept = value),
             linewidth = .01)+
  labs(x = 'Sample size (n)', 
       y = 'Cosine similarity of\nsecond active\nsubspace direction',
       color = 'C matrix rank')

ggplot(data = data.frame(grid_reduced, value = rowMeans(cos_sim3))) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p, labeller = label_both) +
  geom_hline(data = data.frame(grid_reduced, value = rowMeans(cos_sim_ref)), aes(yintercept = value),
             linewidth = .01)+
  labs(x = 'Sample size (n)', 
       y = 'Cosine similarity of\nthird active\nsubspace direction',
       color = 'C matrix rank')




ggplot(data = data.frame(grid_reduced, value = sqrt((sapply(C_mat, function(x) x[2,2]) - 
                                                       sapply(posterior_means, function(x) x[2,2]))^2))) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p, labeller = label_both) +
  scale_y_log10()+ 
  labs(x = 'Sample size (n)', 
       y = 'Posterior mean of C[2,2] RMSE from truth',
       color = 'C matrix rank')
ggplot(data = data.frame(grid_reduced, value = sqrt((sapply(C_mat, function(x) x[1,1]) - 
                                                       sapply(posterior_means, function(x) x[1,1]))^2))) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p, labeller = label_both) +
  scale_y_log10()+
  labs(x = 'Sample size (n)', 
       y = 'Posterior mean of C[1,1] RMSE from truth',
       color = 'C matrix rank')


ggplot(data = data.frame(grid_reduced, value = sqrt((sapply(C_mat, function(x) x[1,2]) - 
                                                       sapply(posterior_means, function(x) x[1,2]))^2))) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p) +
  scale_y_log10()+
  labs(x = 'Sample size (n)', 
       y = 'Posterior mean of C[1,2] RMSE from truth',
       color = 'C matrix rank')

ggplot(data = data.frame(grid_reduced, value = ((sapply(posterior_means, function(x) x[1,1])) - 
                                                  sapply(C_mat, function(x) x[1,1]) 
))) +
  geom_boxplot(aes(x = n, y = value * p, group = paste(n, type), color = type)) +
  facet_wrap(~p, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Posterior mean of C[1,1] difference from truth (bias)',
       color = 'C matrix rank')

### eigenvalues
ggplot(data = data.frame(grid_reduced, value = mean_eigen[,1])) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  # compare with true eigenvalue
  geom_hline(data = data.frame(grid_reduced, value = 500/grid_reduced$p * exp(0)), 
             aes(yintercept = value),
             linewidth = .2) + 
  facet_wrap(~p, labeller = label_both) +
  labs(y = 'Posterior mean first eigenvalue', 
       x = 'Sample size (n)', color = 'C matrix rank')

ggplot(data = data.frame(grid_reduced, value = mean_eigen[,2])) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p) + 
  geom_hline(data = data.frame(grid_reduced, value = 500/grid_reduced$p * exp(-1)), 
             aes(yintercept = value),
             linewidth = .2) + 
  facet_wrap(~p, labeller = label_both) +
  labs(y = 'Posterior mean second eigenvalue', 
       x = 'Sample size (n)', color = 'C matrix rank')

ggplot(data = data.frame(grid_reduced, value = mean_eigen[,3])) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p) + 
  geom_hline(data = data.frame(grid_reduced, value = 500/grid_reduced$p * exp(-2)), 
             aes(yintercept = value),
             linewidth = .2) +
  geom_hline(aes(yintercept = 0),
             linewidth = .2) + 
  facet_wrap(~p, labeller = label_both) +
  labs(y = 'Posterior mean third eigenvalue', 
       x = 'Sample size (n)', color = 'C matrix rank')


ggplot(data = data.frame(grid_reduced, value = mean_eigen[,5])) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p) + 
  geom_hline(data = data.frame(grid_reduced, value = 500/grid_reduced$p * exp(-3)), 
             aes(yintercept = value),
             linewidth = .2) +
  geom_hline(aes(yintercept = 0),
             linewidth = .2) 

ggplot(data = data.frame(grid_reduced, value = mean_eigen[,5])) +
  geom_boxplot(aes(x = n, y = value, group = paste(n, type), color = type)) +
  facet_wrap(~p) + 
  geom_hline(data = data.frame(grid_reduced, value = 500/grid_reduced$p * exp(-4)), 
             aes(yintercept = value),
             linewidth = .2) +
  geom_hline(aes(yintercept = 0),
             linewidth = .2) 




