
d = 4
expected_values <- c(.1, .2, .3, .4)
variances <- c(.04, .02, .02, .02)

MC <- 1000
iter <- 500
mean_theta <- matrix(nrow = d-1, ncol = iter)
var_theta <- matrix(nrow = d, ncol = iter)
mean_theta[,1] <- rep(1, d)
var_theta[,1] <- log(rep(.2, d))
for (i in 2:iter) {
  prev_mean <- mean_theta[,i - 1] 
  prev_var <- var_theta[,i - 1] 
  sample_x <- matrix(nrow = d, ncol = MC)
  sample_x_var <- matrix(nrow = d, ncol = MC)
  sample_U <- matrix(nrow = d, ncol = MC)
  for (mc in 1:MC) {
    sample_x[,mc] <- rnorm(n = d, c(1,prev_mean), sd = sqrt(exp(prev_var)))
    sample_U[,mc] <- (sample_x[,mc] - c(1,prev_mean))/exp(prev_var)
    sample_U_var[,mc] <- - 1/2 +  ((sample_x[,mc] - c(1, prev_mean))^2)/(2*exp(prev_var))
    sample_x[,mc] <- exp(sample_x[,mc] )/sum(exp(sample_x[,mc]))
  }
  mu_hat <- apply(sample_x, 1, mean)
  mu_var_hat <- (apply(sample_x, 1, var))
  
  mu_prime_hat <- apply(sapply(1:ncol(sample_x),
                               function(mc) sample_x[,mc] %*% t(sample_U[,mc]), 
                               simplify = 'array' ),
                        c(1,2), 
                        mean)
  mu_prime_var_hat <- apply(sapply(1:ncol(sample_x),
                                   function(mc) sample_x[,mc] %*% t(sample_U_var[,mc]), 
                               simplify = 'array' ),
                        c(1,2), 
                        mean)
  
  mean_theta[,i] <- prev_mean + solve(crossprod(mu_prime_hat[,2:d]),
                                      crossprod(mu_prime_hat[,2:d], expected_values - mu_hat))
  var_theta[,i] <- prev_var +
                            solve(crossprod(mu_prime_var_hat) + diag(.0001, nrow =d),
                                 crossprod(mu_prime_var_hat, (variances) - mu_var_hat))
  var_theta[,i][var_theta[,i] < -10] <- -10
  var_theta[,i][var_theta[,i] > 10] <- 10
}
cbind(mu_hat, expected_values)
cbind(mu_var_hat, variances)
c(1, mean_theta[,i])
exp(var_theta[,i])
exp(mu_var_hat)
var_theta
plot(var_theta[1,])

plot(mean_theta[1,])
plot(mean_theta[2,])
plot(mean_theta[3,])
plot(var_theta[1,])
plot(var_theta[2,])
plot(var_theta[3,])