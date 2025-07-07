sim.conj_prior = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] prior_scale;
  real <lower=0> prior_dof;
  vector[k] y[N];
}
parameters {
  cov_matrix[k] Q1;
}
transformed parameters {
  matrix[k,k] S;
  matrix[k,k] posterior_scale;
  real <lower=0> posterior_dof;
  
  S = rep_matrix(0.0, k, k); // Initialize S to a zero matrix
  for (n in 1:N) {
    S += y[n] * y[n]'; // Accumulate outer products
  }
  
  posterior_scale <- prior_scale + S;
  posterior_dof <- prior_dof + N;
}
model {
  Q1 ~ inv_wishart(posterior_dof, posterior_scale);
}
"