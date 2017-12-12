import numpy as np 
import pystan as ps
import pickle as pk
from matplotlib import pyplot as plt 

# read in field and computer simulation data
DATACOMP = np.genfromtxt('DATACOMP.csv', delimiter = ',', skip_header=1)
DATAFIELD = np.genfromtxt('DATAFIELD.csv', delimiter = ',', skip_header=1)


y = DATAFIELD[:, 0] # observed output
xf = DATAFIELD[:, 1:] # observed input
(n, p) = xf.shape 
eta = DATACOMP[:, 0] # simulation output
xc = DATACOMP[:, 1:(p+1)] # simulation input
tc = DATACOMP[:, (p+1):] # calibration parameters
(m, q) = tc.shape

x_pred = xf # design points for predictions
n_pred = x_pred.shape[0] # number of predictions

# standardization of output y and eta
eta_mu = np.nanmean(eta) # mean value
eta_sd = np.nanstd(eta) # standard deviation
y = (y - eta_mu) / eta_sd
eta = (eta - eta_mu) /eta_sd

# Put design points xf and xc on [0,1]
x = np.concatenate((xf,xc), axis=0)
x_min = np.nanmin(x, axis=0)
x_max = np.nanmax(x, axis=0)
xf = (xf - x_min) / (x_max - x_min)
xc = (xc - x_min) / (x_max - x_min)
x_pred = (x_pred - x_min) / (x_max - x_min)

# Put calibration parameters t on domain [0,1]
tc_min = np.nanmin(tc, axis=0)
tc_max = np.nanmax(tc, axis=0)
tc = (tc - tc_min) / (tc_max - tc_min)

# create data as dict for input to Stan
stan_data = dict(n=n, m=m, n_pred=n_pred, p=p, y=y, q=q, 
	eta=eta, xf=xf, xc=xc, x_pred=x_pred, tc=tc)

# run model in stan
model = ps.StanModel(file="bcWithPred.stan")
fit = model.sampling(data=stan_data,
	iter=500, 
	chains=3)

# plot traceplots and posterior probability histograms
fit.plot()
plt.savefig("trace_hist_plots.png")

# extract predictions, excluding warm-up and 
y_pred = fit.extract()["y_pred"]
y_pred = y_pred * eta_sd + eta_mu  # convert back to original scale

# for loop to visualize predictions at different input x
for i in range(0, p):
	x_value = np.around(DATAFIELD[:, (i+1)], decimals=1) # round x axis for plotting

	plt.figure(figsize=(9, 6))
	plt.boxplot(y_pred, positions=x_value, 
		flierprops=dict(markersize=0.2))
	plt.scatter(x=x_value, y=DATAFIELD[:, 0], 
		s=10, c="#D55E00")
	plt.savefig("plot%d.png" % (i+1)) # save plot as png file
