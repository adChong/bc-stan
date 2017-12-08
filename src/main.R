library(rstan)
source("util.R")

# read in field and computer simulation data
DATACOMP <- read.csv("DATACOMP.csv", header = TRUE)
DATAFIELD <- read.csv ("DATAFIELD.csv", header = TRUE)

# get dimensions of dataset
p <- ncol(DATAFIELD) - 1 # number of input factors
q <- ncol(DATACOMP) - p - 1 # number of calibration parameters
n <- nrow(DATAFIELD) # sample size of observed field data
m <- nrow(DATACOMP) # sample size of computer simulation data

# extract data from DATAFIELD (Table 3) and DATACOMP (Table 4) 
y <- DATAFIELD[,1] # observed output
xf <- DATAFIELD[,2:(1+p)] # observed input
eta <- DATACOMP[,1] # simulation output
xc <- DATACOMP[,2:(1+p)] # simulation input
tc <- DATACOMP[,(2+p):(1+p+q)] # calibration parameters

x_pred <- xf # design points for predictions
n_pred <- nrow(x_pred) # design points for predictions

# standardization of output y and eta
eta_mu <- mean(eta, na.rm = TRUE) # mean value
eta_sd <- sd(eta, na.rm = TRUE) # standard deviation
y <- (y - eta_mu) / eta_sd
eta <- (eta - eta_mu) / eta_sd

# Put design points xf and xc on [0,1]
x <- rbind(as.matrix(xf), as.matrix(xc))
for (i in (1:ncol(x))){
  x_min <- min(x[,i], na.rm = TRUE)
  x_max <- max(x[,i], na.rm = TRUE)
  xf[,i] <- (xf[,i] - x_min) / (x_max - x_min)
  xc[,i] <- (xc[,i] - x_min) / (x_max - x_min)
  x_pred[,i] <- (x_pred[,i] - x_min) / (x_max - x_min)
}

# Put calibration parameters t on domain [0,1]
for (j in (1:ncol(tc))){
  tc_min <- min(tc[,j], na.rm = TRUE)
  tc_max <- max(tc[,j], na.rm = TRUE)
  tc[,j] <- (tc[,j] - tc_min) / (tc_max - tc_min)
}

# create data as list for input to Stan
stan_data <- list(n=n, m=m, n_pred=n_pred, p=p, y=y, q=q, eta=eta, 
                  xf=as.matrix(xf), xc=as.matrix(xc), 
                  x_pred=as.matrix(x_pred), tc=as.matrix(tc))

# set stan to execute multiple Markov chains in parallel
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# run model in stan
# To run without predictive inferenec: 
# comment lines 58-61 and 85
# uncomment lines 63-66 and 86
fit <- stan(file = "bcWithPred.stan", 
            data = stan_data, 
            iter = 500, 
            chains = 3)

#fit <- stan(file = "bcWithoutPred.stan",
#            data = stan_data,
#            iter = 500,
#            chains = 3)

# plot traceplots, excluding warm-up
stan_trace(fit, pars = c("tf", "beta_eta", "beta_delta", 
                         "lambda_eta", "lambda_delta", "lambda_e"))

# summarize results
print(fit, pars = c("tf", "beta_eta", "beta_delta", 
                    "lambda_eta", "lambda_delta", "lambda_e"))

# posterior probability distribution of tf
stan_hist(fit, pars = c("tf"))

# extract predictions, excluding warm-up and 
samples <- rstan::extract(fit)

# get predictive inference y_pred and convert back to original scale
y_pred <- samples$y_pred * eta_sd + eta_mu 
#y_pred <- y.pred(x_pred, samples, stan_data) * eta_sd + eta_mu 


n_samples <- nrow(y_pred)

# for loop to visualize predictions at different input x
for (i in (1:p)) {
  field_data <- data.frame(yf=DATAFIELD[, 1], 
                           xf=signif(DATAFIELD[, (i+1)], 3))
  pred_data <- matrix(data=t(y_pred), 
                      nrow=length(y_pred), ncol = 1)
  plot_data <- data.frame(apply(field_data, 2, rep, 
                                n_samples), pred = pred_data)
  # save plot as png file
  png(paste("plot", i, ".png", sep = ""))
  plt <- ggplot(data = plot_data, aes(y=pred, x=xf, group=xf)) +
    geom_boxplot(outlier.size=0.2) +
    geom_point(data = field_data, aes(x=xf, y=yf), 
               color="#D55E00", size=0.8) 
  print(plt)
  dev.off()
}







