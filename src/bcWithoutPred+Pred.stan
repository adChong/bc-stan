functions {
  vector gp_pred_rng(int n, int N, int n_pred, int p, int q, vector z, 
                     row_vector[] xf, row_vector[] xc, row_vector[] x_pred,
                     row_vector tf, row_vector[] tc,
                     row_vector beta_eta, real lambda_eta,
                     row_vector beta_delta, real lambda_delta,
                     real delta) {
    vector[n_pred] y_pred;
    {
      matrix[N,N] sigma_eta;  //sigma_eta + sigma_delta = sigma_z
      matrix[n,n] sigma_delta;  //sigma_eta + sigma_delta = sigma_z
      matrix[N,N] sigma_z;  //sigma_z = s11
      matrix[N,N] L_K;  //L_K'*L_K = s11
      vector[N] L_K_div_y;
      vector[N] K_div_y;  //K_div_y = s11{-1} y
      matrix[N, n_pred] k_x1_x2;  //k_x1_x2 = s12
      vector[n_pred] y_pred_mu;  //y_pred_mu = E{y(x_pred)|y,.} = s21 s11{-1} y
      matrix[N,n_pred] v_pred;  //v_pred' * v_pred = s21 s11{-1} s12
      matrix[n_pred,n_pred] sigma_eta_pred; //sigma_eta_pred + sigma_delta_pred = s22
      matrix[n_pred,n_pred] sigma_delta_pred;  //sigma_eta_pred + sigma_delta_pred = s22
      matrix[n_pred,n_pred] cov_y_pred;  //cov_y_pred = COV{y(x_pred, tf)|y,.} = s22 - s21 s11{-1} s12
      row_vector[p+q] xt[N];
      row_vector[p+q] xt_pred[n_pred];
      
      // xt = [[xt,tf],[xc,tc]]
      for (i in 1:n) {
        xt[i] = append_col(xf[i],tf);
      }
      for (i in (n+1):N) {
        xt[i] = append_col(xc[i-n],tc[i-n]);
      }
      for (i in 1:n_pred) {
        xt_pred[i] = append_col(x_pred[i],tf);
      }

	  // Evaluate s11
      for (i in 1:(N-1)) {
        sigma_eta[i, i] = 1/lambda_eta + delta;
        for (j in (i+1):N) {
          sigma_eta[i, j] = exp(-dot_self((xt[i] - xt[j]) .* beta_eta))/lambda_eta;
          sigma_eta[j, i] = sigma_eta[i, j];
        }
      }
      sigma_eta[N, N] = 1/lambda_eta + delta;

      // elements of sigma_delta and add observation errors
      for (i in 1:(n-1)) {
        sigma_delta[i, i] = 1/lambda_delta + 1e-9;
        for (j in (i+1):n) {
          sigma_delta[i, j] = exp(-dot_self((xf[i] - xf[j]) .* beta_delta))/lambda_delta;
          sigma_delta[j, i] = sigma_delta[i, j];
        }
      }
      sigma_delta[n, n] = 1/lambda_delta + 1e-9;

      // computation of covariance matrix sigma_z 
      sigma_z = sigma_eta;
      sigma_z[1:n, 1:n] = sigma_eta[1:n, 1:n] + sigma_delta;

      L_K = cholesky_decompose(sigma_z);

	  // Evaluate E{y(x_pred)|y,.} = s21 s11{-1} y
	  // K_div_y = s11{-1} y
      L_K_div_y = mdivide_left_tri_low(L_K, z);
      K_div_y = mdivide_right_tri_low(L_K_div_y', L_K)';

	  // Evaluate s21 s12
      for (i in 1:N) {
        for (j in 1:n_pred) {
          k_x1_x2[i, j] = exp(-dot_self((xt[i] - xt_pred[j]) .* beta_eta)) / lambda_eta;
        }
      }
      for (i in 1:n) {
        for (j in 1:n_pred) {
          k_x1_x2[i, j] = k_x1_x2[i, j] + exp(-dot_self((xt[i][1:p] - xt_pred[j][1:p]) .* beta_delta)) / lambda_delta;
        }
      }
      y_pred_mu = (k_x1_x2' * K_div_y);
	  
	  // v_pred' * v_pred = s21 s11{-1} s12
      v_pred = mdivide_left_tri_low(L_K, k_x1_x2);

	  // Evaluate s22
      for (i in 1:(n_pred-1)) {
        sigma_eta_pred[i, i] = 1/lambda_eta + delta;
        for (j in (i+1):n_pred) {
          sigma_eta_pred[i, j] = exp(-dot_self((xt_pred[i] - xt_pred[j]) .* beta_eta))/lambda_eta;
          sigma_eta_pred[j, i] = sigma_eta_pred[i, j];
        }
      }
      sigma_eta_pred[n_pred, n_pred] = 1/lambda_eta + delta;

      for (i in 1:(n_pred-1)) {
        sigma_delta_pred[i, i] = 1/lambda_delta + 1e-9;
        for (j in (i+1):n_pred) {
          sigma_delta_pred[i, j] = exp(-dot_self((xt_pred[i][1:p] - xt_pred[j][1:p]) .* beta_delta))/lambda_delta;
          sigma_delta_pred[j, i] = sigma_delta_pred[i, j];
        }
      }
      sigma_delta_pred[n_pred, n_pred] = 1/lambda_delta + 1e-9;

      cov_y_pred = sigma_eta_pred + sigma_delta_pred;
	  
	  // Evaluate COV{y(x_pred, tf)|y,.} = s22 - s21 s11{-1} s12
      cov_y_pred = cov_y_pred - v_pred' * v_pred + diag_matrix(rep_vector(1e-9, n_pred));
      y_pred = multi_normal_rng(y_pred_mu, cov_y_pred);
    }
    return y_pred;
  }
}
data {
  int<lower=0> n; // number of field data
  int<lower=0> m; // number of computer simulation
  int<lower=1> n_pred; // number of predictions
  int<lower=0> p; // number of observable inputs x
  int<lower=0> q; // number of calibration parameters t
  vector[n] y; // field observations
  vector[m] eta; // output of computer simulations
  row_vector[p] xf[n]; // observable inputs corresponding to y
  // (xc, tc): design points corresponding to eta
  row_vector[p] xc[m]; 
  row_vector[q] tc[m]; 
  row_vector[p] x_pred[n_pred]; 
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
  real<lower=0> delta;
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
  matrix[N, N] sigma_eta; // simulator covariance
  matrix[n, n] sigma_delta; // bias term covariance
  matrix[N, N] sigma_z; // covariance matrix
  matrix[N, N] L; // cholesky decomposition of covariance matrix 
  row_vector[p+q] xt[N];
  row_vector[p+q] xt_pred[n_pred];

  // xt = [[xt,tf],[xc,tc]]
  for (i in 1:n) {
    xt[i] = append_col(xf[i],tf);
  }
  for (i in (n+1):N) {
    xt[i] = append_col(xc[i-n],tc[i-n]);
  }
  for (i in 1:n_pred) {
    xt_pred[i] = append_col(x_pred[i],tf);
  }

  // elements of sigma_eta
  for (i in 1:(N-1)) {
	  sigma_eta[i, i] = 1/lambda_eta + delta;
    for (j in (i+1):N) {
	    sigma_eta[i, j] = exp(-dot_self((xt[i] - xt[j]) .* beta_eta))/lambda_eta;
      sigma_eta[j, i] = sigma_eta[i, j];
    }
  }
  sigma_eta[N, N] = 1/lambda_eta + delta;

  // elements of sigma_delta and add observation errors
  for (i in 1:(n-1)) {
	  sigma_delta[i, i] = 1/lambda_delta + 1e-9;
    for (j in (i+1):n) {
	    sigma_delta[i, j] = exp(-dot_self((xf[i] - xf[j]) .* beta_delta))/lambda_delta;
      sigma_delta[j, i] = sigma_delta[i, j];
    }
  }
  sigma_delta[n, n] = 1/lambda_delta + 1e-9;

  // computation of covariance matrix sigma_z 
  sigma_z = sigma_eta;
  sigma_z[1:n, 1:n] = sigma_eta[1:n, 1:n] + sigma_delta;

  // Specify priors here
  rho_eta ~ beta(1.0, 0.3);
  rho_delta ~ beta(1.0, 0.3);
  lambda_eta ~ gamma(10, 10); // gamma (shape, rate)
  lambda_delta ~ gamma(10, 0.3);
  delta ~ normal(0, 1);

  L = cholesky_decompose(sigma_z); // cholesky decomposition
  z ~ multi_normal_cholesky(mu, L);
}

generated quantities {
  vector[n_pred] f_pred = gp_pred_rng(n, N, n_pred, p, q, z, xf, xc, x_pred, tf, tc, beta_eta, lambda_eta, beta_delta, lambda_delta, delta);
  vector[n_pred] y_pred;
  for (i in 1:n_pred) {
    y_pred[i] = normal_rng(f_pred[i], delta);
  }
} 


