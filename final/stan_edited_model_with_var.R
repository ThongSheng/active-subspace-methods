

sim.ss_dirichlet_wishart = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] R;
  vector[N] y;
  vector[N] mu0;
  matrix[N,k] locs;
  vector[k] alpha;
  real <lower=0> prior_cor_dof;
  real <lower=0> prior_gamma_a;
  real <lower=0> prior_gamma_b;
  real <lower=0> prior_sigma2_a;
  real <lower=0> prior_sigma2_b;
  real <lower=0> prior_tau_a;
  real <lower=0> prior_tau_b;
}
parameters {
  cov_matrix[k] Q1;
  simplex[k] xi;
  real <lower=0> K;
  real <lower=0> sigma2;
  real <lower=0> tau;
}
transformed parameters {
  corr_matrix[k] Rho; 
  cov_matrix[k] Sigma;
  cov_matrix[N] Sigma_gp;
  vector<lower=0>[k] delta1;
// Rho is the correlation matrix prior, start with a Q1 ~ IW() and its transformed into
// a correlation matrix with D1*Q1*D1, wehre D1<-diag(delta1), is done with for loops

  for (i in 1:k) delta1[i] <- 1/sqrt(Q1[i,i]);
  for (n in 1:k) {
    for (m in 1:n) {
      Rho[m,n] <- delta1[m] * delta1[n] * Q1[m,n]; 
    }
  }

  for (n in 1:k) {
    for (m in (n+1):k) {
      Rho[m,n] <- Rho[n,m];
    }
  } 

// compute covariance matrix as: Sigma = K * D*Q*D, where D = diag(delta) 
  for (n in 1:k) {
    for (m in 1:n) {
      Sigma[m,n] <- K * sqrt(xi[m]) * sqrt(xi[n]) * Rho[m,n]; 
    }
  }
  for (n in 1:k) {
    for (m in (n+1):k) {
      Sigma[m,n] <- Sigma[n,m];
    }
  }
  
  for (i in 1:N) {
    for (m in 1:i) {
      Sigma_gp[m,i] <- sigma2 * exp(- (locs[i,]- locs[m,]) * Sigma * (locs[i,]- locs[m,])' /2 ); 
    }
    Sigma_gp[i,i] <- Sigma_gp[i,i] + tau;
  }

  for (n in 1:N) {
    for (m in (n+1):N) {
      Sigma_gp[m,n] <- Sigma_gp[n,m];
    }
  } 
}
model {
  Q1 ~ wishart(prior_cor_dof, R);
  xi ~ dirichlet(alpha);
  K ~ gamma(prior_gamma_a, prior_gamma_b);
  sigma2 ~ inv_gamma(prior_sigma2_a, prior_sigma2_b);
  tau ~ inv_gamma(prior_tau_a, prior_tau_b);
  y ~ multi_normal(mu0, Sigma_gp);
}
"

