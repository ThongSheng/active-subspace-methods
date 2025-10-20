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
  
  # Calculate Frobenius Norm between true C and posterior mean of Sigma
  posterior_mean_Sigma <- pred_C
  frobenius <- sqrt(sum((C - posterior_mean_Sigma)^2))
  
  # Calculate cosine similarity of first eigenvector
  C_eigen <- eigen(C)$vectors[,1]
  Sigma_eigen <- eigen(posterior_mean_Sigma)$vector[,1]
  cos_sim_C_Sigma <- abs(sum(C_eigen * Sigma_eigen) / (sqrt(sum(C_eigen^2)) * sqrt(sum(Sigma_eigen^2))))
  
  # Calculate first eigenvalue
  first_eigenvalue <- eigen(posterior_mean_Sigma)$values[1]
  
  # Store results in a list
  nonstan_results_list[[i]] <- data.frame(
    d = grid$p[file_id],
    n = grid$n[file_id],
    seed = grid$seed[file_id],
    func = grid$func[file_id],
    method = grid$prior_choice[file_id],
    time_used = as.double(time_used),
    frobenius = frobenius,
    cos_sim_C_Sigma = cos_sim_C_Sigma,
    first_eigenvalue = first_eigenvalue
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
  
  # Calculate Frobenius Norm between true C and posterior mean of Sigma
  posterior_mean_Sigma <- apply(extract_vals$Sigma, c(2, 3), mean)
  frobenius <- sqrt(sum((C - posterior_mean_Sigma)^2))
  
  # Calculate cosine similarity of first eigenvector
  C_eigen <- eigen(C)$vectors[,1]
  Sigma_eigen <- eigen(posterior_mean_Sigma)$vector[,1]
  cos_sim_C_Sigma <- abs(sum(C_eigen * Sigma_eigen) / (sqrt(sum(C_eigen^2)) * sqrt(sum(Sigma_eigen^2))))
  
  # Calculate first eigenvalue
  first_eigenvalue <- eigen(posterior_mean_Sigma)$values[1]
  
  # Store results in a list
  stan_results_list[[i]] <- data.frame(
    d = grid$p[file_id],
    n = grid$n[file_id],
    seed = grid$seed[file_id],
    func = grid$func[file_id],
    method = grid$prior_choice[file_id],
    time_used = as.double(time_used),
    frobenius = frobenius,
    cos_sim_C_Sigma = cos_sim_C_Sigma,
    first_eigenvalue = first_eigenvalue
  )
}

# Combine the list into a single data frame
stan_results_df <- do.call(rbind, stan_results_list)

# --- Combine both STAN and NONSTAN ---
results_df <- rbind(nonstan_results_df, stan_results_df)

# Plot 1: Computation time
(comp_time <- ggplot(data = results_df) +
    geom_boxplot(aes(x = factor(n), y = time_used, color = method)) +
    facet_grid(d ~ func, labeller = label_both) +
    labs(x = 'Sample size (n)', 
         y = 'Computation time (minutes)') +
    theme_bw())
ggsave("Desktop/time.png", plot = comp_time, width = 8, height = 5, units = "in", dpi = 300)

# Plot 1.5: Log of Computation Time
(log_time <- ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = log(time_used), color = method)) +
  facet_grid(d ~ func, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Log of computation time (minutes)') +
  theme_bw())
ggsave("Desktop/log_time.png", plot = log_time, width = 8, height = 5, units = "in", dpi = 300)

# Plot 2: Frobenius Norm between C and posterior mean of Sigma
(frob_outlier <- ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = frobenius, color = method)) +
  facet_grid(d ~ func, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Frobenius Norm') +
  theme_bw())
ggsave("Desktop/frob_outlier.png", plot = frob_outlier, width = 8, height = 5, units = "in", dpi = 300)

# Plot 2.5: Frobenius Norm without outlier
results_df_clean <- results_df[results_df$frobenius <= 1000,]
(frob <- ggplot(data = results_df_clean) +
    geom_boxplot(aes(x = factor(n), y = frobenius, color = method)) +
    facet_grid(d ~ func, labeller = label_both) +
    labs(x = 'Sample size (n)', 
         y = 'Frobenius Norm') +
    theme_bw())
ggsave("Desktop/frob.png", plot = frob, width = 8, height = 5, units = "in", dpi = 300)

# Plot 3: Cosine Similarity between first eigenvectors
line_data <- data.frame(
  d = c(2, 10, 20),
  yintercept = c(sqrt(2/pi)/sqrt(2), sqrt(2/pi)/sqrt(10), sqrt(2/pi)/sqrt(20))
)

(cosine <- ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = cos_sim_C_Sigma, color = method)) +
  geom_hline(data = line_data, aes(yintercept = yintercept), color = "red", linetype = "dashed") +
  facet_grid(d ~ func, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Cosine Similarity of First Eigenvectors') +
  theme_bw())
ggsave("Desktop/cosine.png", plot = cosine, width = 8, height = 5, units = "in", dpi = 300)

# Plot 4: Posterior mean first eigenvalue
(first_eigen <- ggplot(data = results_df_clean) +
  geom_boxplot(aes(x = factor(n), y = first_eigenvalue, color = method)) +
  facet_grid(d ~ func, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'First Eigenvalue') +
  theme_bw())
ggsave("Desktop/first_eigen.png", plot = first_eigen, width = 8, height = 5, units = "in", dpi = 300)

