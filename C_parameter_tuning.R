library(CholWishart)
library(reshape2)
library(mvtnorm)
library(tictoc)
library(dplyr)
library(stringr)

###############################
#      CBASS (varying n)      #
###############################
set.seed(123)
p <- 2
n_list <- seq(20, 300, 40)
C <- matrix(1/45*c(120, 50, 50, 21), nrow=p)
iters <- 50
CBASS_n <- list()
CBASS_n_time <- list()

for (iter_val in 1:iters) {
  for (n_val in n_list) {
    tic()
    Xunif <- matrix(runif(n_val*p), nrow=n_val, ncol=p)
    
    # create covariance matrix with squared exponential kernel and anisotropy matrix C
    unscaled_dist1 <- matrix(Xunif[,1], nrow = n_val, ncol = n_val) - matrix(Xunif[,1], nrow = n_val, ncol = n_val, byrow = T)
    unscaled_dist2 <- matrix(Xunif[,2], nrow = n_val, ncol = n_val) - matrix(Xunif[,2], nrow = n_val, ncol = n_val, byrow = T)
    dist_sq_mat <- unscaled_dist1^2 *C[1,1] + 
      unscaled_dist2^2 *C[2,2] + 
      unscaled_dist2 * unscaled_dist1 *C[1,2] +
      unscaled_dist2 * unscaled_dist1 *C[2,1]
    Cov_mat_true <- exp(-dist_sq_mat/2)
    diag(Cov_mat_true) <- diag(Cov_mat_true) + .00001 # add jitter to ensure Cov_mat is invertible
    
    # sample data from multivariate normal (assuming the squared-exponential is the "correct" model)
    y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat_true))
    
    # create bass model
    mod_bass <- bass(Xunif, y, verbose=TRUE)
    CBass <- C_bass(mod_bass)
    timing <- toc(quiet = TRUE)
    
    elapsed <- timing$toc - timing$tic
    current_time_record <- data.frame(iter = iter_val, n_val = n_val, time_taken = elapsed)
    CBASS_n_time <- rbind(CBASS_n_time, current_time_record)
    
    temp <- melt(CBass)
    names(temp) <- c("row", "column", "value")
    temp$n <- n_val
    temp$iteration <- iter_val
    CBASS_n <- rbind(CBASS_n, temp)
  }
}


###############################
#   Cconjugate (varying df)   #
###############################
set.seed(123)
p <- 2
n <- 300
df_list <- seq(p+2,300,20)
C <- matrix(1/45*c(120, 50, 50, 21), nrow=p)
iters <- 50
Cconj_df <- list()
Cconj_df_time <- list()

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

for (iter_val in 1:iters) {
  for (df_val in df_list) {
    # Prior parameters
    tic()
    df_prior <- df_val
    Scale_prior <- diag(p) * df_prior
    
    # Generate data and calculate gradient
    Xunif <- matrix(runif(n * p), nrow = n, ncol = p)
    dist_tensor <- array(dim = c(n, n, p))
    for (j in 1:p) {
      dist_tensor[,,j] <- matrix(Xunif[,j], nrow = n, ncol = n) - matrix(Xunif[,j], nrow = n, ncol = n, byrow = T)
    }
    dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
    dist_tensor_mat_reduced_crossprod <- lapply(1:n, 
                                                function(i) lapply(1:n, function(j) {
                                                  tcrossprod(dist_tensor[i,j,])
                                                }))
    Cov_mat <- compute_grad_cov_mat_from_C(C, dist_tensor_mat, n)
    chol_Cov_mat <- chol(Cov_mat)
    y <- t(chol_Cov_mat) %*% rnorm(n * p)
    grad <- matrix(y, nrow = n, ncol = p, byrow = TRUE)
    S <- crossprod(grad)
    
    # Posterior parameters
    Scale_posterior <- Scale_prior + S
    df_posterior <- df_prior + n
    
    # Posterior sampling
    posterior_samples <- rInvWishart(10000, df_posterior, Scale_posterior)
    Cconj <- apply(posterior_samples, 1:2, mean)
    timing <- toc(quiet = TRUE)
    
    # Record time
    elapsed <- timing$toc - timing$tic
    current_time_record <- data.frame(iter = iter_val, df_val = df_val, time_taken = elapsed)
    Cconj_df_time <- rbind(Cconj_df_time, current_time_record)
    
    # Store results
    temp <- melt(Cconj)
    names(temp) <- c("row", "column", "value")
    temp$df <- df_val
    temp$iteration <- iter_val
    Cconj_df <- rbind(Cconj_df, temp)
  }
}


###############################
#    Cconjugate (varying n)   #
###############################
set.seed(123)
p <- 2
n_list <- seq(20, 300, 40)
C <- matrix(1/45*c(120, 50, 50, 21), nrow=p)
iters <- 100
Cconj_n <- list()
Cconj_n_time <- list()

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

for (iter_val in 1:iters) {
  for (n_val in n_list) {
    # Prior parameters
    tic()
    df_prior <- p + 10
    Scale_prior <- diag(p) * df_prior
    
    # Generate data and calculate gradient
    Xunif <- matrix(runif(n_val * p), nrow = n_val, ncol = p)
    dist_tensor <- array(dim = c(n_val, n_val, p))
    for (j in 1:p) {
      dist_tensor[,,j] <- matrix(Xunif[,j], nrow = n_val, ncol = n_val) - matrix(Xunif[,j], nrow = n_val, ncol = n_val, byrow = T)
    }
    dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n_val^2, ncol = p)
    dist_tensor_mat_reduced_crossprod <- lapply(1:n_val, 
                                                function(i) lapply(1:n_val, function(j) {
                                                  tcrossprod(dist_tensor[i,j,])
                                                }))
    Cov_mat <- compute_grad_cov_mat_from_C(C, dist_tensor_mat, n_val)
    chol_Cov_mat <- chol(Cov_mat)
    y <- t(chol_Cov_mat) %*% rnorm(n_val * p)
    grad <- matrix(y, nrow = n_val, ncol = p, byrow = TRUE)
    S <- crossprod(grad)
    
    # Posterior parameters
    Scale_posterior <- Scale_prior + S
    df_posterior <- df_prior + n_val
    
    # Posterior sampling
    posterior_samples <- rInvWishart(10000, df_posterior, Scale_posterior)
    Cconj <- apply(posterior_samples, 1:2, mean)
    timing <- toc(quiet = TRUE)
    
    # Record time
    elapsed <- timing$toc - timing$tic
    current_time_record <- data.frame(iter = iter_val, n_val = n_val, time_taken = elapsed)
    Cconj_n_time <- rbind(Cconj_n_time, current_time_record)
    
    # Store results
    temp <- melt(Cconj)
    names(temp) <- c("row", "column", "value")
    temp$n <- n_val
    temp$iteration <- iter_val
    Cconj_n <- rbind(Cconj_n, temp)
  }
}


########################
#    CGP (varying n)   #
########################
set.seed(123)
p <- 2
C <- matrix(1/45*c(120, 50, 50, 21), nrow=p)
n_list <- seq(20, 300, 40)
lower_lbfgsb <- c(log(0.01), -1, log(0.01))
upper_lbfgsb <- c(log(10), 1, log(10))
iters <- 50
CGP_n <- list()
CGP_n_time <- list()

for (iter_val in 1:iters) {
  for (n_val in n_list) {
    tic()
    Xunif <- matrix(ncol = p, nrow = n_val, runif(p*n_val))
    
    # make matrices of distances in the direction of each dimension
    unscaled_dist1 <- matrix(Xunif[,1], nrow = n_val, ncol = n_val) - matrix(Xunif[,1], nrow = n_val, ncol = n_val, byrow = T)
    unscaled_dist2 <- matrix(Xunif[,2], nrow = n_val, ncol = n_val) - matrix(Xunif[,2], nrow = n_val, ncol = n_val, byrow = T)
    
    # create covariance matrix with squared exponential kernel
    # and anisotropy matrix C
    dist_sq_mat <- unscaled_dist1^2 *C[1,1] + 
      unscaled_dist2^2 *C[2,2] + 
      unscaled_dist2 * unscaled_dist1 *C[1,2] +
      unscaled_dist2 * unscaled_dist1 *C[2,1]
    
    Cov_mat_true <- exp(-dist_sq_mat/2)
    diag(Cov_mat_true) <- diag(Cov_mat_true) + .00001 # add jitter to ensure Cov_mat is invertible
    
    # sample data from multivariate normal
    # we are assuming the squared-exponential is the "correct" model
    y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat_true))
    
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
      C_mat_use <- matrix(nrow = 2, ncol =2, 
                          c(exp(par[1]), 
                            par[2]*sqrt(exp(par[1])*exp(par[3])), 
                            par[2]*sqrt(exp(par[1])*exp(par[3])), 
                            exp(par[3])))
      ll_val <- likelihood_C(C_mat_use)
      -ll_val
    }
    
    # numerically optimize likelihood, including on the log scale for variance parameters
    likelihood_optim <- optim(c(log(3), (0), log(1)), 
                              likelihood_function, 
                              method = "L-BFGS-B", 
                              lower = lower_lbfgsb, 
                              upper = upper_lbfgsb)
    likelihood_est <- matrix(nrow = p, ncol =p, 
                             c(exp(likelihood_optim$par[1]), 
                               likelihood_optim$par[2]*sqrt(exp(likelihood_optim$par[1]))*sqrt(exp(likelihood_optim$par[3])), 
                               likelihood_optim$par[2]*sqrt(exp(likelihood_optim$par[1]))*sqrt(exp(likelihood_optim$par[3])), 
                               exp(likelihood_optim$par[3])))
    timing <- toc(quiet = TRUE)
    
    elapsed <- timing$toc - timing$tic
    current_time_record <- data.frame(iter = iter_val, n_val = n_val, time_taken = elapsed)
    CGP_n_time[[length(CGP_n_time) + 1]] <- current_time_record
    
    temp <- melt(likelihood_est)
    names(temp) <- c("row", "column", "value")
    temp$n <- n_val
    temp$iter <- iter_val
    CGP_n[[length(CGP_n) + 1]] <- temp
  }
}

CGP_n <- do.call(rbind, CGP_n)
CGP_n_time <- do.call(rbind, CGP_n_time)


############################
#    CGP (varying scale)   #
############################
set.seed(123)
p <- 2
C <- matrix(1/45*c(120, 50, 50, 21), nrow=p)
scale_list <- list("identity" = diag(p), 
                   "Ctrue" = C,
                   "Ctrue/10" = C/10, 
                   "Ctrue/5" = C/5,
                   "Ctrue/2" = C/2,
                   "Ctrue*2" = C*2, 
                   "Ctrue*5" = C*5,
                   "Ctrue*10" = C*10)
n <- 300
lower_lbfgsb <- c(log(0.01), -1, log(0.01))
upper_lbfgsb <- c(log(10), 1, log(10))
iters <- 50
CGP_scale <- list()
CGP_scale_time <- list()

for (iter_val in 1:iters) {
  for (scale_val in names(scale_list)) {
    tic()
    Xunif <- matrix(ncol = p, nrow = n, runif(p*n))
    
    # make matrices of distances in the direction of each dimension
    unscaled_dist1 <- matrix(Xunif[,1], nrow = n, ncol = n) - matrix(Xunif[,1], nrow = n, ncol = n, byrow = T)
    unscaled_dist2 <- matrix(Xunif[,2], nrow = n, ncol = n) - matrix(Xunif[,2], nrow = n, ncol = n, byrow = T)
    
    scale_use = scale_list[[scale_val]]
    
    # create covariance matrix with squared exponential kernel
    # and anisotropy matrix C
    dist_sq_mat <- unscaled_dist1^2 * scale_use[1,1] + 
      unscaled_dist2^2 * scale_use[2,2] + 
      unscaled_dist2 * unscaled_dist1 * scale_use[1,2] +
      unscaled_dist2 * unscaled_dist1 * scale_use[2,1]
    
    Cov_mat_true <- exp(- dist_sq_mat/2)
    diag(Cov_mat_true) <- diag(Cov_mat_true) + .00001 # add jitter to ensure Cov_mat is invertible
    
    # sample data from multivariate normal
    # we are assuming the squared-exponential is the "correct" model
    y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat_true))
    
    # implement the likelihood function: data | C, which is multivariate normal
    likelihood_C <- function(scale_use) {
      # for d > 2, we will need to generalize; will need efficient matrix implementation
      dist_sq_mat <- unscaled_dist1^2 * scale_use[1,1] + 
        unscaled_dist2^2 * scale_use[2,2] + 
        unscaled_dist2 * unscaled_dist1 * scale_use[1,2] +
        unscaled_dist2 * unscaled_dist1 * scale_use[2,1]
      Cov_mat <- exp(-dist_sq_mat/2)
      diag(Cov_mat) <- diag(Cov_mat) + .00001
      mvtnorm::dmvnorm(y, sigma = Cov_mat, log = T)
    }
    
    likelihood_function <- function(par) {
      C_mat_use <- matrix(nrow = 2, ncol =2, 
                          c(exp(par[1]), 
                            par[2]*sqrt(exp(par[1])*exp(par[3])), 
                            par[2]*sqrt(exp(par[1])*exp(par[3])), 
                            exp(par[3])))
      ll_val <- likelihood_C(C_mat_use)
      -ll_val
    }
    
    # numerically optimize likelihood, including on the log scale for variance parameters
    likelihood_optim <- optim(c(log(3), (0), log(1)), 
                              likelihood_function, 
                              method = "L-BFGS-B", 
                              lower = lower_lbfgsb, 
                              upper = upper_lbfgsb)
    likelihood_est <- matrix(nrow = p, ncol =p, 
                             c(exp(likelihood_optim$par[1]), 
                               likelihood_optim$par[2]*sqrt(exp(likelihood_optim$par[1]))*sqrt(exp(likelihood_optim$par[3])), 
                               likelihood_optim$par[2]*sqrt(exp(likelihood_optim$par[1]))*sqrt(exp(likelihood_optim$par[3])), 
                               exp(likelihood_optim$par[3])))
    
    # Normalize results
    operator <- str_extract(scale_val, "[*/]")
    number <- as.numeric(str_extract(scale_val, "\\d+$"))
    if (!is.na(operator) && !is.na(number)) {
      if (operator == "/") {
        likelihood_est <- likelihood_est * number
      } else if (operator == "*") {
        likelihood_est <- likelihood_est / number
      }
    }
    
    timing <- toc(quiet = TRUE)
    
    elapsed <- timing$toc - timing$tic
    current_time_record <- data.frame(iter = iter_val, scale_val = scale_val, time_taken = elapsed)
    CGP_scale_time[[length(CGP_scale_time) + 1]] <- current_time_record
    
    temp <- melt(likelihood_est)
    names(temp) <- c("row", "column", "value")
    temp$scale <- scale_val
    temp$iter <- iter_val
    CGP_scale[[length(CGP_scale) + 1]] <- temp
  }
}

CGP_scale <- do.call(rbind, CGP_scale)
CGP_scale_time <- do.call(rbind, CGP_scale_time)


#############################
#    CGP_grad (varying n)   #
#############################
set.seed(123)
p <- 2
C <- matrix(1/45*c(120, 50, 50, 21), nrow=p)
lower_lbfgsb <- c(log(0.01), -1, log(0.01))
upper_lbfgsb <- c(log(10), 1, log(10))
n_list <- seq(20, 180, 40) # reduced sample size for faster computation
iters <- 50
CGP_grad_n <- list()
CGP_grad_n_time <- list()

compute_grad_cov_mat_from_C <- function(C, dist_tensor_mat, n, diagonal_add = .00001) {
  test_distances <- matrix(rowSums((dist_tensor_mat %*% C) * dist_tensor_mat), n, n)
  exp_mat <- exp(- test_distances/2)
  NA_mat <- matrix(NA, p, p)
  
  all_mats <- lapply(1:n, 
                     function(i) {lapply(1:n, function(j) {
                       if (j < i) {
                         return(NA_mat) # do not compute lower triangle, saves computation time
                       } else {
                         mat <- - ((C %*% dist_tensor_mat_reduced_crossprod[[i]][[j]] %*% C - C)*
                                     exp_mat[i,j])
                         if (i == j) {
                           diag(mat) <- diag(mat) + diagonal_add
                         }
                         return(mat)
                       }
                     })
                     })
  # combine all matrices into one pn times pn matrix
  cov_mat <- do.call(rbind, lapply(all_mats, function(x) do.call(cbind, x)))
}

# implement the likelihood function: data | C, which is multivariate normal
likelihood_C_grad <- function(C, dist_tensor_mat, n,
                              diagonal_add) {
  Cov_mat <- compute_grad_cov_mat_from_C(C, dist_tensor_mat, n,
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

for (iter_val in 1:iters) {
  for (n_val in n_list) {
    tic()
    Xunif <- matrix(ncol = p, nrow = n_val, runif(p*n_val))
    
    # We will reorganize the distances theta - theta'
    dist_tensor <- array(dim = c(n_val, n_val, p))
    for (j in 1:p) {
      dist_tensor[,,j] <- matrix(Xunif[,j], nrow = n_val, ncol = n_val) -
        matrix(Xunif[,j], nrow = n_val, ncol = n_val, byrow = T)
    }
    
    # # reshape from n times n times p; to n^2 times p
    dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n_val^2, ncol = p)
    
    # for each i and j, compute p times p matrix of (theta_i - theta_j)(theta_i - theta_j)^T
    # this is used in evaulating the covariance of the gradient
    dist_tensor_mat_reduced_crossprod <- lapply(1:n_val, 
                                                function(i) lapply(1:n_val, function(j) {
                                                  tcrossprod(dist_tensor[i,j,])
                                                }))
    
    Cov_mat <- compute_grad_cov_mat_from_C(C, dist_tensor_mat, n_val)
    chol_Cov_mat <- chol(Cov_mat)
    y <- t(chol_Cov_mat) %*% rnorm(n_val * p)
    
    # will get messier for larger p
    likelihood_function <- function(par, dist_tensor_mat, n,
                                    diagonal_add) {
      C_mat_use <- matrix(nrow = 2, ncol =2, c(exp(par[1]), par[2]*sqrt(exp(par[1]))*sqrt(exp(par[3])), par[2]*sqrt(exp(par[1]))*sqrt(exp(par[3])), exp(par[3])))
      if (min(eigen(C_mat_use)$values) < 10^-5) {
        return(10^6) #check that we have a valid positive-definite matrix
      }
      ll_val <- likelihood_C_grad(C_mat_use, dist_tensor_mat, n,
                                  diagonal_add)
      -ll_val
    }

    # numerically optimize likelihood, including on the log scale for variance parameters
    likelihood_optim <- optim(c(log(3), (0), log(1)), 
                              likelihood_function,
                              dist_tensor_mat = dist_tensor_mat,
                              n = n_val, method = 'L-BFGS-B',
                              lower = lower_lbfgsb, 
                              upper = upper_lbfgsb,
                              diagonal_add = .001)
    likelihood_est <- matrix(nrow = p, ncol =p, 
                           c(exp(likelihood_optim$par[1]), 
                           likelihood_optim$par[2]*sqrt(exp(likelihood_optim$par[1]))*sqrt(exp(likelihood_optim$par[3])), 
                           likelihood_optim$par[2]*sqrt(exp(likelihood_optim$par[1]))*sqrt(exp(likelihood_optim$par[3])), 
                           exp(likelihood_optim$par[3])))
    timing <- toc(quiet = TRUE)
    
    elapsed <- timing$toc - timing$tic
    current_time_record <- data.frame(iter = iter_val, n_val = n_val, time_taken = elapsed)
    CGP_grad_n_time[[length(CGP_grad_n_time) + 1]] <- current_time_record
    
    temp <- melt(likelihood_est)
    names(temp) <- c("row", "column", "value")
    temp$n <- n_val
    temp$iter <- iter_val
    CGP_grad_n[[length(CGP_grad_n) + 1]] <- temp
  }
}

CGP_grad_n <- do.call(rbind, CGP_grad_n)
CGP_grad_n_time <- do.call(rbind, CGP_grad_n_time)


#################################
#    CGP_grad (varying scale)   #
#################################
set.seed(123)
p <- 2
C <- matrix(1/45*c(120, 50, 50, 21), nrow=p)
scale_list <- list("identity" = diag(p), 
                   "Ctrue" = C,
                   "Ctrue/10" = C/10, 
                   "Ctrue/5" = C/5,
                   "Ctrue/2" = C/2,
                   "Ctrue*2" = C*2, 
                   "Ctrue*5" = C*5,
                   "Ctrue*10" = C*10)
n <- 120 # reduced for faster computation
lower_lbfgsb <- c(log(0.01), -1, log(0.01))
upper_lbfgsb <- c(log(10), 1, log(10))
iters <- 30
CGP_grad_scale <- list()
CGP_grad_scale_time <- list()

compute_grad_cov_mat_from_C <- function(C, dist_tensor_mat, n, diagonal_add = .00001) {
  # what we want is the diagonal of 
  # dist_tensor_mat %*% C %*% t(dist_tensor_mat)
  # the rowsums approach computes only the diagonal of this matrix
  test_distances <- matrix(rowSums((dist_tensor_mat %*% C) *
                                     dist_tensor_mat), n, n)
  
  exp_mat <- exp(- test_distances/2)
  NA_mat <- matrix(NA, p, p)
  # for each i,j combination, 
  # Compute p times p matrix
  # -(C * (theta_i - theta_j) *(theta_i - theta_j)^T * C - C) *k(theta_i, theta_j)
  
  all_mats <- lapply(1:n, 
                     function(i) {lapply(1:n, function(j) {
                       if (j < i) {
                         return(NA_mat) # for the lower triangle, do not compute
                         # saves computation time
                       } else {
                         mat <- - ((C %*% dist_tensor_mat_reduced_crossprod[[i]][[j]] %*% C - C)*
                                     exp_mat[i,j])
                         if (i == j) {
                           diag(mat) <- diag(mat) + diagonal_add
                         }
                         return(mat)
                       }
                     })
                     })
  # combine all matrices into one pn times pn matrix
  cov_mat <- do.call(rbind, lapply(all_mats, function(x) do.call(cbind, x)))
}

# implement the likelihood function: data | C, which is multivariate normal
likelihood_C_grad <- function(C, dist_tensor_mat, n,
                              diagonal_add) {
  Cov_mat <- compute_grad_cov_mat_from_C(C, dist_tensor_mat, n,
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

for (iter_val in 1:iters) {
  for (scale_val in names(scale_list)) {
    tic()
    Xunif <- matrix(ncol = p, nrow = n, runif(p*n))
    
    # We will reorganize the distances theta - theta'
    dist_tensor <- array(dim = c(n, n, p))
    for (j in 1:p) {
      dist_tensor[,,j] <- matrix(Xunif[,j], nrow = n, ncol = n) -
        matrix(Xunif[,j], nrow = n, ncol = n, byrow = T)
    }
    
    # # reshape from n times n times p; to n^2 times p
    dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
    
    # for each i and j, compute p times p matrix of (theta_i - theta_j)(theta_i - theta_j)^T
    # this is used in evaulating the covariance of the gradient
    dist_tensor_mat_reduced_crossprod <- lapply(1:n, 
                                                function(i) lapply(1:n, function(j) {
                                                  tcrossprod(dist_tensor[i,j,])
                                                }))
    scale_use = scale_list[[scale_val]]
    Cov_mat <- compute_grad_cov_mat_from_C(scale_use, dist_tensor_mat, n)
    
    
    # sample data from multivariate normal
    # we are assuming the squared-exponential is the "correct" model
    # use cholesky to only use upper triangle
    chol_Cov_mat <- chol(Cov_mat)
    y <- t(chol_Cov_mat) %*% rnorm(n * p)
    
    
    # will get messier for larger p
    likelihood_function <- function(par, dist_tensor_mat, n,
                                    diagonal_add) {
      C_mat_use <- matrix(nrow = 2, ncol =2,
                          c(exp(par[1]),
                            par[2]*sqrt(exp(par[1]))*sqrt(exp(par[3])),
                            par[2]*sqrt(exp(par[1]))*sqrt(exp(par[3])),
                            exp(par[3])))
      if (min(eigen(C_mat_use)$values) < 10^-5) {
        return(10^6) #check that we have a valid positive-definite matrix
      }
      ll_val <- likelihood_C_grad(C_mat_use, dist_tensor_mat, n,
                                  diagonal_add)
      -ll_val
    }
    
    # numerically optimize likelihood, including on the log scale for variance parameters
    likelihood_optim <- optim(c(log(3), (0), log(1)), 
                              likelihood_function,
                              dist_tensor_mat = dist_tensor_mat,
                              n = n, method = 'L-BFGS-B',
                              lower = lower_lbfgsb, 
                              upper = upper_lbfgsb,
                              diagonal_add = .001)
    likelihood_est <- matrix(nrow = p, ncol =p, 
                             c(exp(likelihood_optim$par[1]), 
                               likelihood_optim$par[2]*sqrt(exp(likelihood_optim$par[1]))*sqrt(exp(likelihood_optim$par[3])), 
                               likelihood_optim$par[2]*sqrt(exp(likelihood_optim$par[1]))*sqrt(exp(likelihood_optim$par[3])), 
                               exp(likelihood_optim$par[3])))
    
    # Normalize results
    operator <- str_extract(scale_val, "[*/]")
    number <- as.numeric(str_extract(scale_val, "\\d+$"))
    if (!is.na(operator) && !is.na(number)) {
      if (operator == "/") {
        likelihood_est <- likelihood_est * number
      } else if (operator == "*") {
        likelihood_est <- likelihood_est / number
      }
    }
    
    timing <- toc(quiet = TRUE)
    
    elapsed <- timing$toc - timing$tic
    current_time_record <- data.frame(iter = iter_val, scale_val = scale_val, time_taken = elapsed)
    CGP_grad_scale_time[[length(CGP_grad_scale_time) + 1]] <- current_time_record
    
    temp <- melt(likelihood_est)
    names(temp) <- c("row", "column", "value")
    temp$scale <- scale_val
    temp$iter <- iter_val
    CGP_grad_scale[[length(CGP_grad_scale) + 1]] <- temp
  }
}

CGP_grad_scale <- do.call(rbind, CGP_grad_scale)
CGP_grad_scale_time <- do.call(rbind, CGP_grad_scale_time)


#######################
#   MH (Varying df)   #
#######################
set.seed(123)
p <- 2
n <- 300
prior <- diag(2)
C <- matrix(1/45*c(120, 50, 50, 21), nrow=p)
df_list <- seq(p + 2, 300, 20)
CMH_df <- data.frame(
  row = integer(),
  column = integer(),
  value = numeric(),
  df = numeric(),
  sample_num_in_chain = integer() # To identify the sample number post-burn-in
)
CMH_df_time <- data.frame(df_val = numeric(), time_taken = numeric()) # Initialize for rbind

for (df_val in df_list) {
  tic()
  Xunif <- matrix(ncol = p, nrow = n, runif(p*n))
  unscaled_dist1 <- matrix(Xunif[,1], nrow = n, ncol = n) - matrix(Xunif[,1], nrow = n, ncol = n, byrow = T)
  unscaled_dist2 <- matrix(Xunif[,2], nrow = n, ncol = n) - matrix(Xunif[,2], nrow = n, ncol = n, byrow = T)
  dist_sq_mat <- unscaled_dist1^2 *C[1,1] + 
    unscaled_dist2^2 *C[2,2] + 
    unscaled_dist2 * unscaled_dist1 *C[1,2] +
    unscaled_dist2 * unscaled_dist1 *C[2,1]
  Cov_mat <- exp(- dist_sq_mat/2)
  diag(Cov_mat) <- diag(Cov_mat) + .00001
  
  y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))
  
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
  
  MC <- 10000
  C_MCMC <- array(dim = c(p, p, MC))
  C_MCMC[,,1] <- c(1, 0, 0, 1)
  accept_prob <- rep(0, MC)
  accept <- rep(0, MC)
  
  # prior and likelihood values at initialization
  prior_prev <- CholWishart::dInvWishart(C_MCMC[,,1] / df_val, 
                                         Sigma = prior, df = df_val,
                                         log = T)
  likelihood_value_prev <- likelihood_C(C_MCMC[,,1])
  
  df_prop <- p + 1000 # degrees of freedom for proposal distribution
  
  for (mc in 2:MC) {
    
    # propose next sample
    C_prop <- rWishart(n = 1, df = df_prop, C_MCMC[,,mc - 1])[,,1]/df_prop
    
    # compute log-likelihood/prior for proposed C
    likelihood_value_prop <- likelihood_C(C_prop)
    prior_prop <- CholWishart::dInvWishart(C_prop / df_val, 
                                           Sigma = prior, df = df_val,
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
  
  # remove first samples due to initialization
  burn_in_number <- 1000
  samples_array <- C_MCMC[,,(burn_in_number+1):MC]
  num_retained_samples <- dim(samples_array)[3]
  
  # Store all retained samples
  all_samples_for_current_df_list <- vector("list", num_retained_samples)
  
  for (k in 1:num_retained_samples) {
    one_sample_matrix <- samples_array[,,k]
    melted_sample <- reshape2::melt(one_sample_matrix)
    names(melted_sample) <- c("row", "column", "value")
    melted_sample$df <- df_val
    melted_sample$sample_num_in_chain <- k # k-th sample *after burn-in*
    all_samples_for_current_df_list[[k]] <- melted_sample
  }
  all_samples_for_current_df_df <- do.call(rbind, all_samples_for_current_df_list)
  CMH_df <- rbind(CMH_df, all_samples_for_current_df_df)
  
  # Record time
  timing <- toc(quiet = TRUE)
  elapsed <- timing$toc - timing$tic
  current_time_record <- data.frame(df_val = df_val, time_taken = elapsed)
  CMH_df_time <- rbind(CMH_df_time, current_time_record)
}


######################
#   MH (Varying n)   #
######################
set.seed(123)
p <- 2
C <- matrix(1/45*c(120, 50, 50, 21), nrow=p)
n_list <- seq(20, 300, 40)
prior <- diag(2)
df <- p + 10
CMH_n <- data.frame(
  row = integer(),
  column = integer(),
  value = numeric(),
  n = numeric(),
  sample_num_in_chain = integer() # To identify the sample number post-burn-in
)
CMH_n_time <- data.frame(n_val = numeric(), time_taken = numeric()) # Initialize for rbind

for (n_val in n_list) {
  tic()
  Xunif <- matrix(ncol = p, nrow = n_val, runif(p*n_val))
  unscaled_dist1 <- matrix(Xunif[,1], nrow = n_val, ncol = n_val) - matrix(Xunif[,1], nrow = n_val, ncol = n_val, byrow = T)
  unscaled_dist2 <- matrix(Xunif[,2], nrow = n_val, ncol = n_val) - matrix(Xunif[,2], nrow = n_val, ncol = n_val, byrow = T)
  dist_sq_mat <- unscaled_dist1^2 * C[1,1] + 
    unscaled_dist2^2 * C[2,2] + 
    unscaled_dist2 * unscaled_dist1 * C[1,2] +
    unscaled_dist2 * unscaled_dist1 * C[2,1]
  Cov_mat <- exp(- dist_sq_mat/2)
  diag(Cov_mat) <- diag(Cov_mat) + .00001
  
  y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))
  
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
  
  MC <- 10000
  C_MCMC <- array(dim = c(p, p, MC))
  C_MCMC[,,1] <- c(1, 0, 0, 1)
  accept_prob <- rep(0, MC)
  accept <- rep(0, MC)
  
  # prior and likelihood values at initialization
  prior_prev <- CholWishart::dInvWishart(C_MCMC[,,1] / df, 
                                         Sigma = prior, df = df,
                                         log = T)
  likelihood_value_prev <- likelihood_C(C_MCMC[,,1])
  
  df_prop <- p + 1000 # degrees of freedom for proposal distribution
  
  for (mc in 2:MC) {
    
    # propose next sample
    C_prop <- rWishart(n = 1, df = df_prop, C_MCMC[,,mc - 1])[,,1]/df_prop
    
    # compute log-likelihood/prior for proposed C
    likelihood_value_prop <- likelihood_C(C_prop)
    prior_prop <- CholWishart::dInvWishart(C_prop / df, 
                                           Sigma = prior, df = df,
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
  
  # remove first samples due to initialization
  burn_in_number <- 1000
  samples_array <- C_MCMC[,,(burn_in_number+1):MC]
  num_retained_samples <- dim(samples_array)[3]
  
  # Store all retained samples
  all_samples_for_current_n_list <- vector("list", num_retained_samples)
  
  for (k in 1:num_retained_samples) {
    one_sample_matrix <- samples_array[,,k]
    melted_sample <- reshape2::melt(one_sample_matrix)
    names(melted_sample) <- c("row", "column", "value")
    melted_sample$n <- n_val
    melted_sample$sample_num_in_chain <- k # k-th sample *after burn-in*
    all_samples_for_current_n_list[[k]] <- melted_sample
  }
  all_samples_for_current_n_df <- do.call(rbind, all_samples_for_current_n_list)
  CMH_n <- rbind(CMH_n, all_samples_for_current_n_df)
  
  # Record time
  timing <- toc(quiet = TRUE)
  elapsed <- timing$toc - timing$tic
  current_time_record <- data.frame(n_val = n_val, time_taken = elapsed)
  CMH_n_time <- rbind(CMH_n_time, current_time_record)
}
