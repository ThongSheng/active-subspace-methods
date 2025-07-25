library(mvtnorm)
library(rstan)

# Each prior configuration now only needs:
# - stan_model_name: Path to the compiled STAN model RData file
# - default_params: A list of default prior parameters for this model
# - get_specific_data_params_func: A function to return prior-specific parameters for data_input
# - get_inits_func: A function to return initial values for STAN parameters
# - stan_model_string_var: The name of the model defined in STAN_edited_model.R
prior_configs <- list(
  "dirichlet_wishart" = list(
    stan_model_name = "stan/dirichlet_wishart_model.RData",
    stan_model_string_var = "sim.ss_dirichlet_wishart", 
    default_params = list(
      prior_gamma_a = 7,
      prior_gamma_b = .05
    ),
    get_specific_data_params_func = function(p, config) {
      list(
        alpha = rep(1, p),
        prior_cor_dof = p + 5,
        prior_gamma_a = config$default_params$prior_gamma_a,
        prior_gamma_b = config$default_params$prior_gamma_b
      )
    },
    get_inits_func = function(n_chains, p, config) {
      lapply(1:n_chains, function(x) {
        mat <- matrix(runif(p*p, 0, .5), nrow = p, ncol = p)
        diag(mat) <- 1
        mat <- crossprod(mat)
        list(Q1 = mat, xi = rep(1/p, p), K = runif(1, 0.1, 10))
      })
    }
  ),
  "dirichlet_wishart_reduce" = list(
    stan_model_name = "stan/dirichlet_wishart_reduce_model.RData",
    stan_model_string_var = "sim.ss_dirichlet_wishart_dimreduce", 
    default_params = list(
      prior_gamma_a = 7,
      prior_gamma_b = .05
    ),
    get_specific_data_params_func = function(p, config) {
      list(
        k_reduce = p - 1,
        alpha = rep(1, p),
        prior_cor_dof = p + 50,
        prior_gamma_a = config$default_params$prior_gamma_a,
        prior_gamma_b = config$default_params$prior_gamma_b
      )
    },
    get_inits_func = function(n_chains, p, config) {
      lapply(1:n_chains, function(x) {
        mat <- matrix(runif(p*p, 0, .5), nrow = p, ncol = p)
        diag(mat) <- 1
        mat <- crossprod(mat)
        list(Q1 = mat, xi = rep(1/p, p), K = runif(1, 0.1, 10))
      })
    }
  ),
  "lognormal_inverse_wishart" = list(
    stan_model_name = "stan/lniw_gp_rescale_model.RData",
    stan_model_string_var = "sim.sslniw_gp_rescale",
    default_params = list(
      prior_rescale_mean = 0,
      prior_rescale_var = 1
    ),
    get_specific_data_params_func = function(p, config) {
      list(
        prior_lgn_mean = rep(0, p - 1),
        prior_lgn_var = rep(1, p - 1),
        prior_dof = p + 5,
        prior_rescale_mean = config$default_params$prior_rescale_mean,
        prior_rescale_var = config$default_params$prior_rescale_var
      )
    },
    get_inits_func = function(n_chains, p, config) {
      lapply(1:n_chains, function(x) {
        q1_mat <- matrix(runif(p*p, 0, .5), nrow = p, ncol = p)
        diag(q1_mat) <- 1
        q1_mat <- crossprod(q1_mat)
        list(Q1 = q1_mat, xi = rnorm(p - 1, 0, 0.1), K = rnorm(1, 0, 0.1))
      })
    }
  )
  # To add a new prior:
  # "new_prior_name" = list(
  #   stan_model_name = "stan/new_model.RData",
  #   stan_model_string_var = "sim.new_model_name_variable", # this is the variable name in stan_edited_model.R
  #   default_params = list(...),
  #   get_specific_data_params_func = function(p, config) { list(...) },
  #   get_inits_func = function(n_chains, p, config) { list(...) }
  # )
)

# --- Simulation Grid Definition ---
array_id <- as.integer(Sys.getenv('THISJOBVALUE'))
grid <- expand.grid(
  'p' = c(2, 5, 10, 15, 20),
  'type' = c('full_rank', '2d'),
  'n' = c(25, 75, 125, 175),
  'seed' = 1:20,
  'prior_choice' = names(prior_configs) 
)

# Extract parameters for this specific job array run
p <- grid$p[array_id]
n <- grid$n[array_id]
selected_prior_name <- as.character(grid$prior_choice[array_id])
set.seed(grid$seed[array_id])

# Get the configuration for the selected prior
current_prior_config <- prior_configs[[selected_prior_name]]
if (is.null(current_prior_config)) {
  stop(paste0("Error: Configuration for prior '", selected_prior_name, "' not found in prior_configs."))
}

# --- Data Generation  ---
x_obs <- matrix(ncol = p, nrow = n, runif(p*n))
dist_tensor <- array(dim = c(n, n, p))
for (j in 1:p) {
  dist_tensor[,,j] <- matrix(x_obs[,j], nrow = n, ncol = n) - matrix(x_obs[,j], nrow = n, ncol = n, byrow = T)
}
dist_tensor_mat <- matrix(as.vector(dist_tensor), nrow = n^2, ncol = p)
dist_tensor_mat_reduced <- dist_tensor_mat[upper.tri(matrix(nrow =n, ncol = n), diag = T), ]

compute_cov_mat_from_C <- function(C, dist_tensor_mat_reduced, n, diagonal_add = .00001) {
  test_distances <- rowSums((dist_tensor_mat_reduced %*% C) * dist_tensor_mat_reduced)
  dist_mat_C <- matrix(nrow = n, ncol = n)
  dist_mat_C[upper.tri(dist_mat_C, diag = T)] <- test_distances
  cov_mat <- exp(-dist_mat_C/2)
  diag(cov_mat) <- diag(cov_mat) + diagonal_add
  cov_mat
}

W_random <- eigen(crossprod(matrix(rnorm(p * p), nrow = p, ncol = p)))$vectors
if (grid$type[array_id] == 'full_rank') {
  C <- 500/p * (W_random %*% diag(exp((p:1 -p))) %*% t(W_random))
} else {
  C <- 500/p * (W_random %*% diag(exp(c(0, -1, rep(-Inf, p-2)))) %*% t(W_random))
}
Cov_mat <- compute_cov_mat_from_C(C, dist_tensor_mat_reduced, n)
Cov_mat[lower.tri(Cov_mat)] <- t(Cov_mat)[lower.tri(Cov_mat)]
y <- as.vector(rmvnorm(n = 1, sigma = Cov_mat))


# --- STAN Model Loading ---
stan_model_object <- NULL
stan_model_path <- current_prior_config$stan_model_name
stan_code_file_master <- "stan/stan_edited_model.R"

# Ensure the 'stan' directory exists
if (!dir.exists("stan")) {
  dir.create("stan")
}

if (!file.exists(stan_model_path)) {
  if (file.exists(stan_code_file_master)) {
    source(stan_code_file_master)
  } else {
    stop(paste0("Master STAN model code file not found: ", stan_code_file_master))
  }
  
  # Get the specific model string using its variable name
  model_string_var_name <- current_prior_config$stan_model_string_var
  if (!exists(model_string_var_name)) {
    stop(paste0("STAN model string variable '", model_string_var_name, "' not found in ", stan_code_file_master))
  }
  model_string <- get(model_string_var_name)
  
  stan_model_object <- stan_model(model_code = model_string)
  save(stan_model_object, file = stan_model_path)
} else {
  load(stan_model_path)
}


# Construct data_input dynamically
get_data_input <- function(N, k, R, y, mu0, locs, diag_add, specific_params) {
  base_data <- list(
    N = N,
    k = k,
    R = R,
    y = y,
    mu0 = mu0,
    locs = locs,
    diag_add = diag_add
  )
  # Combine base parameters with prior-specific parameters
  c(base_data, specific_params)
}

# Generate prior-specific data parameters using the function from prior_configs
specific_data_params <- current_prior_config$get_specific_data_params_func(p, current_prior_config)
data_input <- get_data_input(n, 
                             p, 
                             diag(x = 1, nrow = p), 
                             y, 
                             rep(0, n), 
                             x_obs, 
                             .00001, 
                             specific_data_params)


# MCMC initializations
n_chains <- 4
it <- 1500
w <- 500

initial_values <- current_prior_config$get_inits_func(n_chains, p, current_prior_config)
sapply(initial_values, function(x) min(eigen(x$Q1)$values))


# Run STAN sampling
a_time <- Sys.time()
out_model <- sampling(object = stan_model_object, data = data_input,
                      pars = c('Sigma'), iter = it, chains = n_chains, warmup = w,
                      init = initial_values,
                      cores = n_chains)
b_time <- Sys.time()
time_used <- b_time - a_time
units(time_used) <- 'mins'

extract_vals <- extract(out_model)
summary_vals <- summary(out_model)

# Save results
output_dir <- '/scratch/negishi/angt/simulation_results/'
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
save_file_name <- paste0(output_dir, selected_prior_name, "_", array_id, '.RData')
save(extract_vals, summary_vals, C, y, time_used, file = save_file_name)
q(save = 'no')
