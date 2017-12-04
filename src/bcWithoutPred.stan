data {
  int<lower=0> n; // number of field data
  int<lower=0> m; // number of computer simulation
  int<lower=0> p; // number of observable inputs x
  int<lower=0> q; // number of calibration parameters t
  vector[n] y; // field observations
  vector[m] eta; // output of computer simulations
  matrix[n, p] xf; // observable inputs corresponding to y
  // (xc, tc): design points corresponding to eta
  matrix[m, p] xc; 
  matrix[m, q] tc; 
}

transformed data {
  int<lower=1> N;
  vector[n+m] z; // z = [y, eta]
  vector[n+m] mu; // mean vector

  N = n + m;
  // set mean vector to zero
  for (i in 1:N) {
    mu[i] = 0;
  }
  z = append_row(y, eta); 
}

parameters {
  // tf: calibration parameters
  // rho_eta: reparameterization of beta_eta
  // rho_delta: reparameterization of beta_delta
  // lambda_eta: precision parameter for eta
  // lambda_delta: precision parameter for bias term
  // lambda_e: precision parameter of observation error
  row_vector<lower=0,upper=1>[q] tf; 
  row_vector<lower=0,upper=1>[p+q] rho_eta;
  row_vector<lower=0,upper=1>[p] rho_delta;
  real<lower=0> lambda_eta; 
  real<lower=0> lambda_delta; 
  real<lower=0> lambda_e; 
}

transformed parameters {
  // beta_delta: correlation parameter for bias term
  // beta_e: correlation parameter of observation error
  row_vector[p+q] beta_eta;
  row_vector[p] beta_delta;
  beta_eta = -4.0 * log(rho_eta);
  beta_delta = -4.0 * log(rho_delta);
}

model {
  // declare variables
  matrix[N, (p+q)] xt;
  matrix[N, N] sigma_eta; // simulator covariance
  matrix[n, n] sigma_delta; // bias term covariance
  matrix[N, N] sigma_z; // covariance matrix
  matrix[N, N] L; // cholesky decomposition of covariance matrix 
  row_vector[p] temp_delta;
  row_vector[p+q] temp_eta;

  // xt = [[xt,tf],[xc,tc]]
  xt[1:n, 1:p] = xf;
  xt[(n+1):N, 1:p] = xc;
  xt[1:n, (p+1):(p+q)] = rep_matrix(tf, n);
  xt[(n+1):N, (p+1):(p+q)] = tc;

  // diagonal elements of sigma_eta
  sigma_eta = diag_matrix(rep_vector((1 / lambda_eta), N));
 
  // off-diagonal elements of sigma_eta
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      temp_eta = xt[i] - xt[j];
      sigma_eta[i, j] = beta_eta .* temp_eta * temp_eta';
      sigma_eta[i, j] = exp(-sigma_eta[i, j]) / lambda_eta;
      sigma_eta[j, i] = sigma_eta[i, j];
    }
  }

  // diagonal elements of sigma_delta
  sigma_delta = diag_matrix(rep_vector((1 / lambda_delta), n));

  // off-diagonal elements of sigma_delta
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      temp_delta = xf[i] - xf[j];
      sigma_delta[i, j] = beta_delta .* temp_delta * temp_delta';
      sigma_delta[i, j] = exp(-sigma_delta[i, j]) / lambda_delta;
      sigma_delta[j, i] = sigma_delta[i, j];
    }   
  }

  // computation of covariance matrix sigma_z 
  sigma_z = sigma_eta;
  sigma_z[1:n, 1:n] = sigma_eta[1:n, 1:n] + sigma_delta;

  // add observation errors
  for (i in 1:n) {
    sigma_z[i, i] = sigma_z[i, i] + (1.0 / lambda_e);
  }  

  // Specify priors here
  rho_eta[1:(p+q)] ~ beta(1.0, 0.3);
  rho_delta[1:p] ~ beta(1.0, 0.3);
  lambda_eta ~ gamma(10, 10); // gamma (shape, rate)
  lambda_delta ~ gamma(10, 0.3); 
  lambda_e ~ gamma(10, 0.03); 

  L = cholesky_decompose(sigma_z); // cholesky decomposition
  z ~ multi_normal_cholesky(mu, L);
}

