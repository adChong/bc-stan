library(fields)
# duplicate functions ###########################################################

rep.row <- function(x, n){
  # Function to create a matrix by repeating the vector n times by row.
  #
  # Args:
  #   x: a vector. 
  #   n: number of times to repeat the vector.
  #
  # Returns:
  #   output: matrix containing vector x repeated n times by row
  
  output <- matrix(rep(x, each=n), nrow=n)
  return(output)
}

# Error Calculation ###########################################################

cvrmse <- function(sim,obs){
  # Function to compute the Coefficient of Variation of the 
  # Root Mean Squared Error (CVRMSE) between two vectors 
  # according to ASHRAE guideline 14.
  #
  # Args:
  #   sim: Vector containing the predicted output.
  #   obs: Vector containing the measured or observed output.
  #
  # Returns:
  #   cvrmse: CVRMSE according to ASHRAE guideline 14.
  n <- length(sim)
  cvrmse <- 100*(sqrt(sum((sim-obs)^2,na.rm=TRUE)/(n-1))/mean(obs,na.rm=TRUE))
  return(cvrmse)
}

nmbe <- function(sim,obs){
  # Functin that computes the Normalized Mean Bias Error (NMBE) between 
  # two vectors according to ASHRAE guideline 14.
  #
  # Args:
  #   sim: Vector containing the predicted output.
  #   obs: Vector containing the measured or observed output.
  #
  # Returns:
  #   nbme: NMBE according to ASHRAE guideline 14.
  n <- length(sim)
  nmbe <- 100*(sum(sim-obs,na.rm=TRUE)/((n-1)*mean(obs,na.rm=TRUE)))
  return(nmbe)
}


# Distance Calculation ###########################################################
calc.distance <- function(x, alpha=2){
  # Function to create a distance structure
  #
  # Args:
  #   x: a numeric matrix. 
  #   alpha: power of the Minkowski distance (generally assume alpha = 2).
  #
  # Returns:
  #   dist_lower: vector containing off-diagonal lower triangle entries 
  #               of the distance matrix for x
  
  dist <- rdist(x)^alpha
  dist_lower <- dist[lower.tri(dist)]
  return(dist_lower)
}

gen.dist <- function(x){
  n <- nrow(x) # number of data
  # create distance structure containing pariwise distance computation for 
  # each column (default value of alpha = 2)
  d <- apply(x, 2, calc.distance, alpha=2)
  return(list(n=n,d=as.matrix(d)))
}

gasp.cov <- function(x, beta, lam){
  # Function to generates the Gaussian process (GP) covariance matrix sigma
  #
  # Args:
  #   x: a numeric matrix (n x p)  
  #   beta: GP correlation hyperparameter
  #   lam: GP precision hyperparameter
  #
  # Returns:
  #   GP covariance matrix sigma (n x n)
  #   sigma = 1/lam exp{-sum_{k=1:p} beta(k)*(x(k)-x'(k))^2} 
  #   where
  #   n: number of rows in x
  #   p: number of columns in x
  
  n <- nrow(x) # number of data
  # create distance structure containing pariwise distance computation for 
  # each column (default value of alpha = 2)
  d <- as.matrix(apply(x, 2, calc.distance)) 
  
  sigma <- matrix(data=0, nrow=n, ncol=n)
  sigma[lower.tri(sigma, diag=FALSE)] <- exp(-d %*% beta) / lam
  sigma <- sigma + t(sigma)
  diag(sigma) <- 1 / lam
  return(sigma)
}

# Prediction Calculation ###########################################################

y.pred <- function(x_pred, samples, params){
  # Function to compute the predictions y at new (x_pred) values conditional 
  # on the observed/measured data
  #
  # Args:
  #   x_pred: design points to compute the predictions
  #   samples: list containing the following
  #     n: number of field data
  #     m: number of computer simulation data
  #     q: number of observable inputs x
  #     p: number of calibration parameters t
  #     y: field observations / measured data
  #     eta: output of computer simulations
  #     xf: observable inputs corresponding to y
  #     (xc, tc): design points corresponding to eta
  #   params: list containing posterior samples of the following parameters 
  #     tf: calibration parameters      
  #     beta_eta: correlation parameter for eta
  #     beta_delta: correlation parameter for bias term
  #     lambda_eta: precision parameter for eta
  #     lambda_delta: precision parameter for bias term
  #     lambda_e: precision parameter for observation errors
  #
  # Returns:
  #     y_pred: predictions for x_pred sampled from mvrnorm where
  #     c(y, eta, y_pred) ~ N(0, sigma_y_pred)
  #     sigma_y_pred = sigma_eta + ...
  #     | sigma_b    0    sigma_bb* | +  | sigma_e   0     0  |
  #     |    0       0        0     |    |   0       0     0  |
  #     | sigma_b*b  0    sigma_b*  |    |   0       0     0  |
  
  cat('starting run \n')
  q <- params$q
  p <- params$p
  n <- params$n
  m <- params$m
  xf <- params$xf
  xc <- params$xc
  tc <- params$tc
  y <- params$y
  eta <- params$eta
  nm <- n + m
  pq <- p + q
  z <- array(data=NA, dim=c(nm, pq))
  z[1:n, 1:p] <- as.matrix(xf)
  z[(n+1):nm, 1:p] <- as.matrix(xc)
  z[(n+1):nm, (p+1):pq] <- as.matrix(tc)
  
  tf <- samples$tf
  beta_eta <- samples$beta_eta
  beta_delta <- samples$beta_delta
  lambda_eta <- samples$lambda_eta
  lambda_delta <- samples$lambda_delta
  lambda_e <- samples$lambda_e
  
  # Number of realizations at which to evaluate E{eta(x_pred,tf)|y,.}
  numrealize <- length(lambda_e)
  
  # Number of predictions
  if(p==1){
    npred <- length(x_pred)
    xx <- 1
  }else{
    npred <- nrow(x_pred)
    xx <- ncol(x_pred)
  }
  
  # Index of y values
  idxy <- c((1:n),(nm+1):(nm+npred))
  
  # storage
  y_pred <- array(data=NA, dim=c(numrealize, npred)) # eta + delta
  
  
  for (I in (1:numrealize)){
    cat('simulation ',I,' of ',numrealize,'\n')
    
    z[1:n, (p+1):(p+q)] <- as.matrix(rep.row(tf[I,], n))
    xtpred <- cbind(as.matrix(x_pred),
                    rep.row(tf[I,], npred))
    
    # Compile design points and prediction points from computer experiments
    Z <- rbind(z,as.matrix(xtpred))
    # Evaluate sigma_eta
    sigma_eta <- gasp.cov(Z, beta_eta[I,], lambda_eta[I])
    # Evaluate sigma_b
    sigma_b <- gasp.cov(xf, beta_delta[I,], lambda_delta[I])
    
    # Compile design points xf and prediction points x_pred for y
    X <- if (p==1) c(xf, x_pred) else rbind(xf, x_pred)
    
    # Evaluate sigma_B
    sigma_B <- gasp.cov(X, beta_delta[I,], lambda_delta[I])
    
    # Evaluate sigma_e
    sigma_e <- diag(n) * (1 / lambda_e[I])
    
    # Compile sigma_y_pred
    sigma_y_pred <- sigma_eta
    sigma_y_pred[idxy,idxy] <- sigma_y_pred[idxy,idxy] + sigma_B
    sigma_y_pred[1:n,1:n] <- sigma_y_pred[1:n,1:n] + sigma_e
    
    
    # solve covariance 
    # Evaluate s11
    s11 <- sigma_y_pred[1:nm, 1:nm]
    # Evaluate s21
    s21 <- sigma_y_pred[(nm+1):nrow(sigma_y_pred), 1:nm]
    # Evaluate s22
    s22 <- sigma_y_pred[(nm+1):nrow(sigma_y_pred), (nm+1):ncol(sigma_y_pred)]
    # Evaluate s12
    s12 <- sigma_y_pred[1:nm, (nm+1):ncol(sigma_y_pred)]
    
    # Evaluate E{y(x_pred)|y,.} = s21 s11{-1} y
    s11y <- solve(s11, c(y, eta))
    E <- s21 %*% s11y
    
    # Evaluate COV{y(x_pred, tf)|y,.} = s22 - s21 s11{-1} s12
    s11s12 <- solve(s11, s12)
    s21s11s12 <- s21 %*% s11s12
    C <- s22 - s21s11s12
    
    # Generate a MVN sample
    y_pred[I,] <- mvrnorm(1, E, C)
    
  }
  cat('run ended \n')
  return(y_pred)
}

eta.pred <- function(x_pred, samples, params){
  # Function to compute the predictions eta (y - delta) at new (x_pred) values 
  # conditional on the observed/measured data
  #
  # Args:
  #   x_pred: design points to compute the predictions
  #   samples: list containing the following
  #     n: number of field data
  #     m: number of computer simulation data
  #     q: number of observable inputs x
  #     p: number of calibration parameters t
  #     y: field observations / measured data
  #     eta: output of computer simulations
  #     xf: observable inputs corresponding to y
  #     (xc, tc): design points corresponding to eta
  #   params: list containing posterior samples of the following parameters 
  #     tf: calibration parameters      
  #     beta_eta: correlation parameter for eta
  #     beta_delta: correlation parameter for bias term
  #     lambda_eta: precision parameter for eta
  #     lambda_delta: precision parameter for bias term
  #     lambda_e: precision parameter for observation errors
  #
  # Returns:
  #     y_pred: predictions for x_pred sampled from mvrnorm where
  #     c(y, eta, y_pred) ~ N(0, sigma_eta_pred)
  #     sigma_eta_pred = sigma_eta + ...
  #     | sigma_b   0   sigma_bb* | +  | sigma_e   0     0  |
  #     |     0     0       0     |    |   0       0     0  |
  #     |     0     0       0     |    |   0       0     0  |

  cat('starting run \n')
  q <- params$q
  p <- params$p
  n <- params$n
  m <- params$m
  xf <- params$xf
  xc <- params$xc
  tc <- params$tc
  y <- params$y
  eta <- params$eta
  nm <- n+m
  pq <- p+q
  z <- array(data=NA, dim=c(nm, pq))
  z[1:n, 1:p] <- as.matrix(xf)
  z[(n+1):nm, 1:p] <- as.matrix(xc)
  z[(n+1):nm, (p+1):pq] <- as.matrix(tc)
  
  tf <- samples$tf
  beta_eta <- samples$beta_eta
  beta_delta <- samples$beta_delta
  lambda_eta <- samples$lambda_eta
  lambda_delta <- samples$lambda_delta
  lambda_e <- samples$lambda_e
  
  # Number of realizations at which to evaluate E{eta(x_pred,theta)|y,.}
  num_realize <- length(lambda_e)
  
  # Number of predictions
  if(p==1){
    npred <- length(x_pred)
    xx <- 1
  }else{
    npred <- nrow(x_pred)
    xx <- ncol(x_pred)
  }
  
  # Index of y values
  idxy <- c((1:n),(nm+1):(nm+npred))
  
  # storage
  eta_pred <- array(data=NA, dim=c(numrealize,npred)) # eta 
  y_pred <- array(data=NA, dim=c(numrealize,npred)) # eta + delta
  
  
  for (I in (1:numrealize)){
    cat('simulation ',I,' of ',numrealize,'\n')
    
    z[1:n,(p+1):(p+q)] <- as.matrix(rep.row(tf[I,],n))
    xtpred <- cbind(as.matrix(x_pred),
                    rep.row(tf[I,],npred))
    
    # Compile design points and prediction points from computer experiments
    Z <- rbind(z,as.matrix(xtpred))
    # Evaluate sigma_eta
    sigma_eta <- gasp.cov(Z,beta_eta[I,],lambda_eta[I])
    # Evaluate sigma_b
    sigma_b <- gasp.cov(xf,beta_delta[I,],lambda_delta[I])
    
    # Compile dsign points xf and prediction points x_pred for y
    X <- if (p==1) c(xf, x_pred) else rbind(xf, x_pred)
    
    # Evaluate sigma_B
    sigma_B <- gasp.cov(X, beta_delta[I,], lambda_delta[I])
    
    # Evaluate sigma_e
    sigma_e <- diag(n) * (1/lambda_e[I])
    
    # Compile covariance matrix sigma_eta_pred
    sigma_eta_pred <- sigma_eta
    sigma_eta_pred[1:n, 1:n] <- sigma_eta[1:n, 1:n] + sigma_b + sigma_e
    
    # solve covariance 
    # Evaluate s11
    s11 <- sigma_eta_pred[1:nm, 1:nm]
    # Evaluate s21
    s21 <- sigma_eta_pred[(nm+1):nrow(sigma_eta_pred), 1:nm]
    # Evaluate s22
    s22 <- sigma_eta_pred[(nm+1):nrow(sigma_eta_pred), 
                          (nm+1):ncol(sigma_eta_pred)]
    # Evaluate s12
    s12 <- sigma_eta_pred[1:nm, (nm+1):ncol(sigma_eta_pred)]
    
    # Evaluate E{y(x_pred)|y,.} = s21 s11{-1} y
    s11y <- solve(s11, c(y,eta))
    E <- s21 %*% s11y
    
    # Evaluate COV{y(x_pred,theta)|y,.} = s22 - s21 s11{-1} s12
    s11s12 <- solve(s11, s12)
    s21s11s12 <- s21 %*% s11s12
    C <- s22 - s21s11s12
    
    eta_pred[I,] <- mvrnorm(1, E, C)
    
    
  }
  cat('run ended \n')
  return(eta_pred)
  
}



