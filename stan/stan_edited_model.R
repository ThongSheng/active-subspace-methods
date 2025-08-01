# implements separated inverse wishart for cov matrix of MVN
# from Alvarez et al
sim.sslniw = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] R;
  vector[k] y[N];
  vector[k] mu0;
  real <lower=0> prior_dof;
  vector[k] prior_lgn_mean;
  vector[k] prior_lgn_var;
}
parameters {
//  vector[k] mu;
  cov_matrix[k] Q1;
  vector[k] xi;
}
transformed parameters {
  corr_matrix[k] Rho; 
  cov_matrix[k] Sigma;
  vector<lower=0>[k] delta;
  vector<lower=0>[k] delta1;
  real s1;
  real s2;
  real rho;
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

// compute covariance matrix as: Sigma = D*Q*D, where D = diag(delta) 
  for (i in 1:k)  delta[i] <- exp( xi[i] );
  for (n in 1:k) {
    for (m in 1:n) {
      Sigma[m,n] <- delta[m] * delta[n] * Rho[m,n]; 
    }
  }
  for (n in 1:k) {
    for (m in (n+1):k) {
      Sigma[m,n] <- Sigma[n,m];
    }
  }
  s1 <- sqrt( Sigma[1,1]);
  s2 <- sqrt(Sigma[2,2]) ;
  rho <-  Sigma[1,2] /(s1*s2);
}
model {
    Q1 ~ inv_wishart(prior_dof, R);
    for ( i in 1:k) {
        xi[i] ~ normal(prior_lgn_mean[k], prior_lgn_var[k]);
    }
  for (n in 1:N) y[n] ~ multi_normal(mu0, Sigma);
}
"

# implements scaled prior for cov matrix of MVN
# adapted from Alvarez et al
sim.sslniw_rescale = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] R;
  vector[k] y[N];
  vector[k] mu0;
  real <lower=0> prior_dof;
  vector[k-1] prior_lgn_mean;
  vector[k-1] prior_lgn_var;
  real prior_rescale_mean;
  real prior_rescale_var;
}
parameters {
//  vector[k] mu;
  cov_matrix[k] Q1;
  vector[k] xi;
  real K;
}
transformed parameters {
  corr_matrix[k] Rho; 
  cov_matrix[k] Sigma;
  vector<lower=0>[k] delta;
  vector<lower=0>[k] delta1;
  real s1;
  real s2;
  real rho;
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

// compute covariance matrix as: Sigma = D*Q*D, where D = diag(delta) 
  for (i in 1:(k-1))  delta[i] <- exp( xi[i] );
  delta[k] <- exp( 1 - sum(xi));
  
  for (n in 1:k) {
    for (m in 1:n) {
      Sigma[m,n] <- exp(K) * delta[m] * delta[n] * Rho[m,n]; 
    }
  }
  for (n in 1:k) {
    for (m in (n+1):k) {
      Sigma[m,n] <- Sigma[n,m];
    }
  }
  s1 <- sqrt(Sigma[1,1]);
  s2 <- sqrt(Sigma[2,2]) ;
  rho <-  Sigma[1,2] /(s1*s2);
}
model {
    Q1 ~ inv_wishart(prior_dof, R);
    for (i in 1:(k-1)) {
        xi[i] ~ normal(prior_lgn_mean[i], sqrt(prior_lgn_var[i]));
    }
    K ~ normal(prior_rescale_mean, sqrt(prior_rescale_var));
  for (n in 1:N) y[n] ~ multi_normal(mu0, Sigma);
}
"

# implements GP approach without rescaling
# adapted from Alvarez et al
sim.sslniw_gp = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] R;
  vector[N] y;
  vector[N] mu0;
  matrix[N,k] locs;
  vector[k] prior_lgn_mean;
  vector[k] prior_lgn_var;
  real <lower=0> prior_dof;
  real <lower=0> diag_add;
}
parameters {
  cov_matrix[k] Q1;
  vector[k] xi;
}
transformed parameters {
  corr_matrix[k] Rho; 
  cov_matrix[k] Sigma;
  cov_matrix[N] Sigma_gp;
  vector<lower=0>[k] delta;
  vector<lower=0>[k] delta1;
  real s1;
  real s2;
  real rho;
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

// compute covariance matrix as: Sigma = D*Q*D, where D = diag(delta) 
  for (i in 1:k)  delta[i] <- exp( xi[i] );
  for (n in 1:k) {
    for (m in 1:n) {
      Sigma[m,n] <- delta[m] * delta[n] * Rho[m,n]; 
    }
  }
  for (n in 1:k) {
    for (m in (n+1):k) {
      Sigma[m,n] <- Sigma[n,m];
    }
  }
  
  for (i in 1:N) {
    for (m in 1:i) {
      Sigma_gp[m,i] <- exp(- (locs[i,]- locs[m,]) * Sigma * (locs[i,]- locs[m,])' /2 ); 
    }
    Sigma_gp[i,i] <- Sigma_gp[i,i] + diag_add;
  }

  for (n in 1:N) {
    for (m in (n+1):N) {
      Sigma_gp[m,n] <- Sigma_gp[n,m];
    }
  } 
}
model {
    Q1 ~ inv_wishart(prior_dof, R);
    for ( i in 1:k) {
        xi[i] ~ normal(prior_lgn_mean[i], sqrt(prior_lgn_var[i]));
    }
  y ~ multi_normal(mu0, Sigma_gp);
}
"


# implements propoosed approach and prior
# adapted from Alvarez et al
sim.sslniw_gp_rescale = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] R;
  vector[N] y;
  vector[N] mu0;
  matrix[N,k] locs;
  vector[k-1] prior_lgn_mean;
  vector[k-1] prior_lgn_var;
  real <lower=0> prior_dof;
  real prior_rescale_mean;
  real prior_rescale_var;
  real <lower=0> diag_add;
}
parameters {
  cov_matrix[k] Q1;
  vector[k-1] xi;
  real K;
}
transformed parameters {
  corr_matrix[k] Rho; 
  cov_matrix[k] Sigma;
  cov_matrix[N] Sigma_gp;
  vector<lower=0>[k] delta;
  vector<lower=0>[k] delta1;
  real s1;
  real s2;
  real rho;
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
  for (i in 1:(k-1))  delta[i] <- exp( xi[i] );
  delta[k] <- exp( 1 - sum(xi));
  for (n in 1:k) {
    for (m in 1:n) {
      Sigma[m,n] <- exp(K) * delta[m] * delta[n] * Rho[m,n]; 
    }
  }
  for (n in 1:k) {
    for (m in (n+1):k) {
      Sigma[m,n] <- Sigma[n,m];
    }
  }
  s1 <- sqrt(Sigma[1,1]);
  s2 <- sqrt(Sigma[2,2]) ;
  rho <-  Sigma[1,2] /(s1*s2);
  
  for (i in 1:N) {
    for (m in 1:i) {
      Sigma_gp[m,i] <- exp(- (locs[i,]- locs[m,]) * Sigma * (locs[i,]- locs[m,])' /2 ); 
    }
    Sigma_gp[i,i] <- Sigma_gp[i,i] + diag_add;
  }

  for (n in 1:N) {
    for (m in (n+1):N) {
      Sigma_gp[m,n] <- Sigma_gp[n,m];
    }
  } 
}
model {
    Q1 ~ inv_wishart(prior_dof, R);
    for ( i in 1:(k-1)) {
        xi[i] ~ normal(prior_lgn_mean[i], sqrt(prior_lgn_var[i]));
    }
    K ~ normal(prior_rescale_mean, sqrt(prior_rescale_var));
  y ~ multi_normal(mu0, Sigma_gp);
}
"


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
  real <lower=0> diag_add;
}
parameters {
  cov_matrix[k] Q1;
  simplex[k] xi;
  real <lower=0> K;
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
      Sigma_gp[m,i] <- exp(- (locs[i,]- locs[m,]) * Sigma * (locs[i,]- locs[m,])' /2 ); 
    }
    Sigma_gp[i,i] <- Sigma_gp[i,i] + diag_add;
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
  y ~ multi_normal(mu0, Sigma_gp);
}
"


sim.ss_dirichlet_wishart_dimreduce = "
data {
  int <lower=0> N;
  int <lower=0> k;
  int <lower=0> k_reduce;
  matrix[k,k] R;
  vector[N] y;
  vector[N] mu0;
  matrix[N,k] locs;
  vector[k] alpha;
  real <lower=0> prior_cor_dof;
  real prior_gamma_a;
  real prior_gamma_b;
  real <lower=0> diag_add;
}
parameters {
  cov_matrix[k] Q1;
  simplex[k] xi;
  real <lower=0> K;
}
transformed parameters {
  matrix[k,k] Rho_prelim = rep_matrix(0, k, k); 
  matrix[k,k] Rho; 
  matrix[k, k] Sigma;
  cov_matrix[N] Sigma_gp;
  vector<lower=0>[k] delta1;
// Rho is the correlation matrix prior, start with a Q1 ~ IW() and its transformed into
// a correlation matrix with D1*Q1*D1, wehre D1<-diag(delta1), is done with for loops
  matrix[k, k] U = eigenvectors_sym(Q1);
  vector[k] Lambda = eigenvalues_sym(Q1);
  
  for (n in 1:k_reduce) {
    Rho_prelim = Rho_prelim + Lambda[n] * col(U,n) * col(U,n)'; 
  }

  for (i in 1:k) delta1[i] <- 1/sqrt(Rho_prelim[i,i]);
  for (n in 1:k) {
    for (m in 1:n) {
      Rho[m,n] <- delta1[m] * delta1[n] * Rho_prelim[m,n]; 
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
      Sigma_gp[m,i] <- exp(- (locs[i,]- locs[m,]) * Sigma * (locs[i,]- locs[m,])' /2 ); 
    }
    Sigma_gp[i,i] <- Sigma_gp[i,i] + diag_add;
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
  y ~ multi_normal(mu0, Sigma_gp);
}
"

