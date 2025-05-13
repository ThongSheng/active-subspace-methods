library(CholWishart)
library(ggplot2)
library(reshape2)

# Define f(x) and its gradient
f <- function(x) x[1]^2 + x[1]*x[2] + x[2]^3/9
grad_f <- function(x) c(2*x[1] + x[2], x[1] + (1/3)*x[2]^2)

# Sim Study Parameters
N <- 1000
p <- 2

# Ctrue
Ctrue <- matrix(1/45*c(120, 50, 50, 21), nrow=p)

# Cbass
set.seed(123)
Xunif <- matrix(runif(N * p), nrow=N, ncol=p)
Yf <- apply(Xunif, 1, f)
mod_bass <- BASS::bass(Xunif, Yf, verbose=FALSE)
Cbass <- C_bass(mod_bass)

# Cbayesian with v=4, n=1000, and 2x2 Identity as prior of scale matrix
Scale_prior <- diag(2)
v_prior <- p + 2
prior_samples <- rInvWishart(N, v_prior, Scale_prior)

X <- matrix(runif(N * p), nrow=N, ncol=p)
grad <- t(apply(X, 1, grad_f))
S <- crossprod(grad)

Scale_posterior <- Scale_prior + S
v_posterior <- v_prior + N

posterior_samples <- rInvWishart(N, v_posterior, Scale_posterior)
Cbayesian <- apply(posterior_samples, 1:2, mean)
Cbayesian

## Testing out Cbayesian with different parameters
# Initialize parameters
Cbayesian_v <- Cbayesian_n <- Cbayesian_prior <- data.frame()
v_list <- seq(p+2,300,20) # df shouldn't be less than 4
n_list <- seq(10, 100, 20) # sample sizes from 10 to 100 with increments of 10
prior_list <- list("identity" = diag(p), "identity_2" = 2*diag(p), "Ctrue" = Ctrue, "Cbass" = Cbass) # priors for scale matrix
iters <- 200



# Changing degrees of freedom (fix n = 1000 and Scale_prior = 2x2 identity)
set.seed(123)
Cbayesian_v <- data.frame()
for (iter in 1:iters) {
  for (v in v_list) {
    # Prior parameters
    v_prior <- v
    Scale_prior <- diag(p) * v_prior
    
    # Generate data and calculate gradient
    X <- matrix(runif(N * p), nrow = N, ncol = p)
    grad <- t(apply(X, 1, grad_f))
    S <- crossprod(grad)
    
    # Posterior parameters
    Scale_posterior <- Scale_prior + S
    v_posterior <- v_prior + N
    
    # Posterior sampling and compute Cbayesian
    posterior_samples <- rInvWishart(N, v_posterior, Scale_posterior)
    Cbayesian <- apply(posterior_samples, 1:2, mean)
    
    # Store results
    temp_df <- melt(Cbayesian)
    names(temp_df) <- c("row", "column", "value")
    temp_df$v <- v
    temp_df$iteration <- iter
    Cbayesian_v <- rbind(Cbayesian_v, temp_df)
  }
}

# Create violin plot
entry_row <- 1; entry_col <- 2
ggplot(subset(Cbayesian_v, row == entry_row & column == entry_col), 
       aes(x = factor(v), y = value, fill = factor(v))) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_hline(yintercept = Ctrue[entry_row, entry_col],
             col = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Degrees of Freedom Prior", 
       y = "Entry value", 
       title = sprintf("Distribution of [%s,%s] Matrix Entry Across Iterations", entry_row, entry_col),
       subtitle = sprintf("True value: %.2f (red dashed line)", Ctrue[entry_row, entry_col])) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d()



# Changing sample sizes (fix v = p+2 and Scale_prior = 2x2 identity)
set.seed(123)
Cbayesian_n <- data.frame()
for (iter in 1:iters) {
  for (n in n_list) {
    # Prior parameters
    v_prior <- p+2
    Scale_prior <- diag(p) * v_prior
    
    # Generate data and calculate gradient
    X <- matrix(runif(n * p), nrow = n, ncol = p)
    grad <- t(apply(X, 1, grad_f))
    S <- crossprod(grad)
    
    # Posterior parameters
    Scale_posterior <- Scale_prior + S
    v_posterior <- v_prior + n
    
    # Posterior sampling and compute Cbayesian
    posterior_samples <- rInvWishart(n, v_posterior, Scale_posterior)
    Cbayesian <- apply(posterior_samples, 1:2, mean)
    
    # Store results
    temp_df <- melt(Cbayesian)
    names(temp_df) <- c("row", "column", "value")
    temp_df$n <- n
    temp_df$iteration <- iter
    Cbayesian_n <- rbind(Cbayesian_n, temp_df)
  }
}

# Create violin plot
entry_row <- 1; entry_col <- 2
ggplot(subset(Cbayesian_n, row == entry_row & column == entry_col), 
       aes(x = factor(n), y = value, fill = factor(n))) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_hline(yintercept = Ctrue[entry_row, entry_col],
             col = "red", linetype = "dashed", linewidth = 1) +
  labs(x = "Sample Size Prior", 
       y = "Entry value", 
       title = sprintf("Distribution of [%s,%s] Matrix Entry Across Iterations", entry_row, entry_col),
       subtitle = sprintf("True value: %.2f (red dashed line)", Ctrue[entry_row, entry_col])) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d()



# Changing scale matrix prior (fix n = 1000 and v = p+2)
set.seed(123)
Cbayesian_prior <- data.frame()
for (iter in 1:iters) {
  for (prior in names(prior_list)) {
    # Prior parameters
    v_prior <- p+200
    Scale_prior <- prior_list[[prior]] * v_prior
    
    # Generate data and calculate gradient
    X <- matrix(runif(N * p), nrow = N, ncol = p)
    grad <- t(apply(X, 1, grad_f))
    S <- crossprod(grad)
    
    # Posterior parameters
    Scale_posterior <- Scale_prior + S
    v_posterior <- v_prior + N
    
    # Posterior sampling and compute Cbayesian
    posterior_samples <- rInvWishart(N, v_posterior, Scale_posterior)
    Cbayesian <- apply(posterior_samples, 1:2, mean)
    
    # Store results
    temp_df <- melt(Cbayesian)
    names(temp_df) <- c("row", "column", "value")
    temp_df$prior <- prior
    temp_df$iteration <- iter
    Cbayesian_prior <- rbind(Cbayesian_prior, temp_df)
  }
}

# Create violin plot
entry_row <- 1; entry_col <- 1
ggplot(subset(Cbayesian_prior, row == entry_row & column == entry_col), 
       aes(x = factor(prior), y = value, fill = factor(prior))) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_hline(yintercept = Ctrue[entry_row, entry_col],
             col = "red", linetype = "dashed", linewidth = 1) +
  
  labs(x = "Scale Matrix Prior", 
       y = "Entry value", 
       title = sprintf("Distribution of [%s,%s] Matrix Entry Across Iterations", entry_row, entry_col),
       subtitle = sprintf("True value: %.2f (red dashed line)", Ctrue[entry_row, entry_col])) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_d()