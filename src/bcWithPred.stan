data {
  int<lower=1> n; // number of field data
  int<lower=1> m; // number of computer simulation
  int<lower=1> n_pred; // number of predictions
  int<lower=1> p; // number of observable inputs x
  int<lower=1> q; // number of calibration parameters t
  vector[n] y; // field observations
  vector[m] eta; // output of computer simulations
  matrix[n, p] xf; // observable inputs corresponding to y
  // (xc, tc): design points corresponding to eta
  matrix[m, p] xc; 
  matrix[m, q] tc; 
  // x_pred: new design points for predictions
  matrix[n_pred, p] x_pred; 
}

transformed data {
  int<lower = 1> N;
  vector[n+m] y_eta;
  vector[n+m+n_pred] mu; // mean vector
  matrix[n+n_pred, p] X; // X=[xf, x_pred]
  
  N = n + m + n_pred;
  // set mean vector to zero
  for (i in 1:N) {
    mu[i] = 0;
  }
  X = append_row(xf, x_pred);
  y_eta = append_row(y, eta);
}

parameters {
  // tf: calibration parameters
  // rho_eta: reparameterization of beta_eta
  // rho_delta: reparameterization of beta_delta
  // lambda_eta: precision parameter for eta
  // lambda_delta: precision parameter for bias term
  // lambda_e: precision parameter of observation error
  // y_pred: predictions
  row_vector<lower=0, upper=1>[q] tf; 
  row_vector<lower=0, upper=1>[p+q] rho_eta; 
  row_vector<lower=0, upper=1>[p] rho_delta; 
  real<lower=0> lambda_eta; 
  real<lower=0> lambda_delta;
  real<lower=0> lambda_e; 
  vector[n_pred] y_pred; 
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
  matrix[N, p+q] xt; 
  matrix[N, N] sigma_eta; // simulator covarinace
  matrix[n+n_pred, n+n_pred] sigma_delta; // bias term covariance
  matrix[n, n] sigma_y; // observation errors
  matrix[N, N] sigma_z; // covariance matrix
  matrix[N, N] L; // cholesky decomposition of covariance matrix 
  vector[N] z; // z = [y, eta, y_pred]
  row_vector[p] temp_delta;
  row_vector[p+q] temp_eta;

  z = append_row(y_eta, y_pred); // z = [y, eta, y_pred]

  // xt = [[xf,tf],[xc,tc],[x_pred,tf]]
  xt[1:n, 1:p] = xf;
  xt[1:n, (p+1):(p+q)] = rep_matrix(tf, n);
  xt[(n+1):(n+m), 1:p] = xc;
  xt[(n+1):(n+m), (p+1):(p+q)] = tc;
  xt[(n+m+1):N, 1:p] = x_pred;
  xt[(n+m+1):N, (p+1):(p+q)] = rep_matrix(tf, n_pred);
  
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
  sigma_delta = diag_matrix(rep_vector((1 / lambda_delta), 
    n+n_pred));
  
  // off-diagonal elements of sigma_delta
  for (i in 1:(n+n_pred-1)) {
    for (j in (i+1):(n+n_pred)) {
      temp_delta = X[i] - X[j];
      sigma_delta[i, j] = beta_delta .* temp_delta * temp_delta';
      sigma_delta[i, j] = exp(-sigma_delta[i, j]) / lambda_delta;
      sigma_delta[j, i] = sigma_delta[i, j];
    }   
  }

  // observation errors sigma_y
  sigma_y = diag_matrix(rep_vector((1 / lambda_e), n));

  // computation of covariance matrix sigma_z 
  sigma_z = sigma_eta;
  sigma_z[1:n, 1:n] = sigma_eta[1:n, 1:n] + 
    sigma_delta[1:n, 1:n] + sigma_y;
  sigma_z[1:n, (n+m+1):N] = sigma_eta[1:n, (n+m+1):N] + 
    sigma_delta[1:n, (n+1):(n+n_pred)];
  sigma_z[(n+m+1):N, 1:n] = sigma_eta[(n+m+1):N, 1:n] + 
    sigma_delta[(n+1):(n+n_pred),1:n];
  sigma_z[(n+m+1):N, (n+m+1):N] = sigma_eta[(n+m+1):N, (n+m+1):N] + 
    sigma_delta[(n+1):(n+n_pred), (n+1):(n+n_pred)];

  // Specify priors here
  rho_eta[1:(p+q)] ~ beta(1.0, 0.3);

  rho_delta[1:p] ~ beta(1.0, 0.3);

  lambda_eta ~ gamma(10, 10); // gamma (shape, rate)

  lambda_delta ~ gamma(10, 0.3); // gamma (shape, rate)
  
  lambda_e ~ gamma(10, 0.03); // gamma (shape, rate)

  L = cholesky_decompose(sigma_z); // cholesky decomposition 
  z ~ multi_normal_cholesky(mu, L);
}
