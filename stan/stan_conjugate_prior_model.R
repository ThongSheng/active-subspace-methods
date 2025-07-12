sim.conj_prior = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] prior_scale;
  real <lower=0> prior_dof;
  vector[k] mu0;
  vector[k] y[N];
}
parameters {
  cov_matrix[k] Q1;
}
model {
  Q1 ~ inv_wishart(prior_dof, prior_scale);
  for (i in 1:N) y[i] ~ multi_normal(mu0, Q1);
}
"