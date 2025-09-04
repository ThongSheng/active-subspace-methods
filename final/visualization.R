# --- Analysis and Visualization ---
library(ggplot2)
library(reshape2)

# --- NON-STAN ----
grid <- expand.grid(
  'p' = c(2, 10, 20),
  'n' = c(20, 100, 150),
  'seed' = 1:3,
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
  
  # Calculate cosine similarity
  C_vec <- as.vector(C)
  Sigma_vec <- as.vector(posterior_mean_Sigma)
  cos_sim_C_Sigma <- abs(sum(C_vec * Sigma_vec) / (sqrt(sum(C_vec^2)) * sqrt(sum(Sigma_vec^2))))
  
  # Calculate RMSE
  rmse <- sqrt(mean((C_vec - Sigma_vec)^2))
  
  # Calculate first eigenvalue
  first_eigenvalue <- eigen(posterior_mean_Sigma)$values[1]
  
  # Store results in a list
  nonstan_results_list[[i]] <- data.frame(
    p = grid$p[file_id],
    n = grid$n[file_id],
    seed = grid$seed[file_id],
    func = grid$func[file_id],
    prior_choice = grid$prior_choice[file_id],
    time_used = as.double(time_used),
    frobenius = frobenius,
    cos_sim_C_Sigma = cos_sim_C_Sigma,
    rmse = rmse,
    first_eigenvalue = first_eigenvalue
  )
}

# Combine the list into a single data frame
nonstan_results_df <- do.call(rbind, nonstan_results_list)


# --- STAN ---
grid <- expand.grid(
  'p' = c(2, 10, 20),
  'n' = c(20, 100, 150),
  'seed' = 1:3,
  'func' = c('determ_full', 'determ_2d', 'GP_full', 'GP_2d'),
  'prior_choice' = c("dirichlet_wishart", "dirichlet_wishart_reduce", "lognormal_inverse_wishart")
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
  
  # Calculate cosine similarity
  C_vec <- as.vector(C)
  Sigma_vec <- as.vector(posterior_mean_Sigma)
  cos_sim_C_Sigma <- abs(sum(C_vec * Sigma_vec) / (sqrt(sum(C_vec^2)) * sqrt(sum(Sigma_vec^2))))
  
  # Calculate RMSE
  rmse <- sqrt(mean((C_vec - Sigma_vec)^2))
  
  # Calculate first eigenvalue
  first_eigenvalue <- eigen(posterior_mean_Sigma)$values[1]
  
  # Store results in a list
  stan_results_list[[i]] <- data.frame(
    p = grid$p[file_id],
    n = grid$n[file_id],
    seed = grid$seed[file_id],
    func = grid$func[file_id],
    prior_choice = grid$prior_choice[file_id],
    time_used = as.double(time_used),
    frobenius = frobenius,
    cos_sim_C_Sigma = cos_sim_C_Sigma,
    rmse = rmse,
    first_eigenvalue = first_eigenvalue
  )
}

# Combine the list into a single data frame
stan_results_df <- do.call(rbind, stan_results_list)


# --- Combine both STAN and NONSTAN ---
results_df <- rbind(nonstan_results_df, stan_results_df)

# Plot 1: Time to compute
ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = time_used, color = prior_choice)) +
  facet_grid(func ~ p, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Time to compute (minutes)', 
       title = 'Computation Time') +
  theme_bw()

# Plot 2: Frobenius Norm between C and posterior mean of Sigma
ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = frobenius, color = prior_choice)) +
  facet_grid(func ~ p, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Frobenius Norm', 
       title = 'Frobenius Norm of Predicted Sigma vs. True C') +
  theme_bw()

# Plot 3: Cosine Similarity between C and posterior mean of Sigma
line_data <- data.frame(
  p = c(2, 10, 20),
  yintercept = c(sqrt(2/pi)/sqrt(2), sqrt(2/pi)/sqrt(10), sqrt(2/pi)/sqrt(20))
)

ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = cos_sim_C_Sigma, color = prior_choice)) +
  geom_hline(data = line_data, aes(yintercept = yintercept), color = "red", linetype = "dashed") +
  facet_grid(func ~ p, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'Cosine Similarity',
       title = 'Cosine Similarity of Predicted Sigma vs. True C') +
  theme(legend.position = 'bottom') +
  theme_bw()

# Plot 4: Posterior mean of predicted value RMSE from truth
ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = rmse, color = prior_choice)) +
  facet_grid(func ~ p, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'RMSE',
       title = 'Posterior Mean of Predicted Value RMSE from Truth') +
  theme_bw()

# Plot 5: Posterior mean first eigenvalue
ggplot(data = results_df) +
  geom_boxplot(aes(x = factor(n), y = first_eigenvalue, color = prior_choice)) +
  facet_grid(func ~ p, labeller = label_both) +
  labs(x = 'Sample size (n)', 
       y = 'First Eigenvalue',
       title = 'Posterior Mean First Eigenvalue') +
  theme_bw()

#results_df <- results_df[results_df$frobenius <= 1e15,] # to remove outliers
