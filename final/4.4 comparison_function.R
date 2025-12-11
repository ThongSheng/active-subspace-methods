# --- Analysis and Visualization ---
library(ggplot2)
library(reshape2)

# --- NON-STAN ----
grid <- expand.grid(
  'p' = c(2, 10, 20),
  'n' = c(20, 100, 150),
  'seed' = 1:30,
  'func' = c('determ_full', 'determ_2d', 'GP_full', 'GP_2d'),
  'prior_choice' = c('bass', 'wycoff', 'mle')
)

output_dir <- '/scratch/negishi/angt/nonstan_results/' # change this to match
files <- list.files(output_dir, pattern = "\\.RData$", full.names = TRUE)
file_ids <- as.integer(gsub(".RData$", "", gsub(".*_", "", files)))

nonstan_results_list <- list()

for (i in 1:length(files)) {
  file_id <- file_ids[i]
  
  # Load the file
  load(files[i])
  
  # Normalize C matrices
  C_norm <- C/sum(diag(C))
  if(sum(diag(pred_C)) == 0) {
    C_est_norm <- NA
  } else {
    C_est_norm <- pred_C / sum(diag(pred_C))
  }
  
  # Calculate Frobenius Norm between normalized C matrices
  frobenius <- sqrt(sum((C_norm - C_est_norm)^2))
  
  # Calculate cosine similarity of first eigenvectors
  C_eigenvector <- eigen(C)$vectors[,1]
  C_est_eigenvector <- eigen(pred_C)$vector[,1]
  cos_sim <- abs(sum(C_eigenvector * C_est_eigenvector))
  
  # Calculate differences in first eigenvalues between normalized matrices
  C_eigenvalue <- ifelse(sum(diag(C)) == 0, NA, eigen(C_norm)$values[1])
  C_est_eigenvalue <- ifelse(sum(diag(pred_C)) == 0, NA, eigen(C_est_norm)$values[1])
  diff_eigenvalue <- abs(C_eigenvalue - C_est_eigenvalue)
  
  # Store results in a list
  nonstan_results_list[[i]] <- data.frame(
    d = grid$p[file_id],
    n = grid$n[file_id],
    seed = grid$seed[file_id],
    func = grid$func[file_id],
    method = grid$prior_choice[file_id],
    time_used = as.double(time_used),
    frobenius = frobenius,
    cos_sim = cos_sim,
    diff_eigenvalue = diff_eigenvalue,
    n_eff = NA
  )
}

# Combine the list into a single data frame
nonstan_results_df <- do.call(rbind, nonstan_results_list)


# --- STAN ---
grid <- expand.grid(
  'p' = c(2, 10, 20),
  'n' = c(20, 100, 150),
  'seed' = 1:30,
  'func' = c('determ_full', 'determ_2d', 'GP_full', 'GP_2d'),
  'prior_choice' = c("dirichlet_wishart", "gp_reduce", "lognormal_inverse_wishart")
)

output_dir <- '/scratch/negishi/angt/stan_results/' # change this to match
files <- list.files(output_dir, pattern = "\\.RData$", full.names = TRUE)
file_ids <- as.integer(gsub(".RData$", "", gsub(".*_", "", files)))

stan_results_list <- list()

for (i in 1:length(files)) {
  file_id <- file_ids[i]
  
  # Load the file
  load(files[i])
  
  # Normalized C matrices
  C_norm <- C/sum(diag(C))
  pred_C <- apply(extract_vals$Sigma, c(2, 3), mean)
  if(sum(diag(pred_C)) == 0) {
    C_est_norm <- NA
  } else {
    C_est_norm <- pred_C / sum(diag(pred_C))
  }
  
  # Calculate Frobenius Norm between normalized C matrices
  frobenius <- sqrt(sum((C_norm - C_est_norm)^2))
  
  # Calculate cosine similarity of first eigenvectors
  C_eigenvector <- eigen(C)$vectors[,1]
  C_est_eigenvector <- eigen(pred_C)$vector[,1]
  cos_sim <- abs(sum(C_eigenvector * C_est_eigenvector))
  
  # Calculate differences in first eigenvalues of normalized C matrices
  C_eigenvalue <- ifelse(sum(diag(C)) == 0, NA, eigen(C_norm)$values[1])
  C_est_eigenvalue <- ifelse(sum(diag(pred_C)) == 0, NA, eigen(C_est_norm)$values[1])
  diff_eigenvalue <- abs(C_eigenvalue - C_est_eigenvalue)
  
  # Additional piece for STAN: n_eff
  n_eff <- mean(head(summary_vals$summary[,9], -1)) # mean of all n_eff
  
  # Store results in a list
  stan_results_list[[i]] <- data.frame(
    d = grid$p[file_id],
    n = grid$n[file_id],
    seed = grid$seed[file_id],
    func = grid$func[file_id],
    method = grid$prior_choice[file_id],
    time_used = as.double(time_used),
    frobenius = frobenius,
    cos_sim = cos_sim,
    diff_eigenvalue = diff_eigenvalue,
    n_eff = n_eff
  )
}

# Combine the list into a single data frame
stan_results_df <- do.call(rbind, stan_results_list)

# --- Combine both STAN and NONSTAN ---
results_df <- rbind(nonstan_results_df, stan_results_df)

# Plot 1: Log of Computation Time
(log_time <- ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = log(time_used), color = method)) +
  facet_grid(d ~ func, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Log of computation time (minutes)') +
  theme_bw() +
  theme(text = element_text(family = 'Arial')))
ggsave("Desktop/log_time.png", plot = log_time, width = 10, height = 5, units = "in", dpi = 300)

# Plot 2: Frobenius Norm between normalized C and C_est
(frob <- ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = frobenius, color = method)) +
  facet_grid(d ~ func, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Frobenius Norm') +
  theme_bw() +
  theme(text = element_text(family = 'Arial')))
ggsave("Desktop/frob.png", plot = frob, width = 10, height = 5, units = "in", dpi = 300)

# Plot 3: Cosine Similarity between first eigenvectors
line_data <- data.frame(
  d = c(2, 10, 20),
  yintercept = c(sqrt(2/pi)/sqrt(2), sqrt(2/pi)/sqrt(10), sqrt(2/pi)/sqrt(20))
)

(cosine <- ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = cos_sim, color = method)) +
  geom_hline(data = line_data, aes(yintercept = yintercept), color = "red", linetype = "dashed") +
  facet_grid(d ~ func, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Cosine Similarity of First Eigenvectors') +
  theme_bw() +
  theme(text = element_text(family = 'Arial')))
ggsave("Desktop/cosine.png", plot = cosine, width = 10, height = 5, units = "in", dpi = 300)

# Plot 4: Differences in first eigenvalue
(first_eigen <- ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = diff_eigenvalue, color = method)) +
  facet_grid(d ~ func, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Abs Difference in First Eigenvalues') +
  theme_bw() +
  theme(text = element_text(family = 'Arial')))
ggsave("Desktop/first_eigen.png", plot = first_eigen, width = 10, height = 5, units = "in", dpi = 300)

# Plot 5: Effective sample size
(ess <- ggplot(data = results_df) +
    geom_boxplot(aes(x = factor(n), y = n_eff, color = method)) +
    facet_grid(d ~ func, labeller = label_both) +
    labs(x = 'Sample size (n)', 
         y = 'Effective sample size (out of 4000)') +
    ylim(0, 4000) +
    theme_bw() +
    theme(text = element_text(family = 'Arial')))
ggsave("Desktop/ess.png", plot = ess, width = 10, height = 5, units = "in", dpi = 300)

