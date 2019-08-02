# bc-stan

## Description
bc-stan provides code for:

1. [Bayesian calibration of Building Energy Models using the Stan modeling language](#Bayesian-calibration-of-energy-models-using-Stan)
2. [Parameter screening with Morris method](#Parameter-screening-using-Morris-method)

## Related Publications
1. Detailed description of this code can be found in [Chong and Menberg (2018)](https://doi.org/10.1016/j.enbuild.2018.06.028)
2. Application case study of Bayesian calibration [Chong et. al (2017)](https://doi.org/10.1016/j.enbuild.2017.08.069)

## Bayesian calibration of energy models using Stan

The stan code for Bayesian calibration is based on [Kennedy and O'Hagan's (2001)](http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00294) Bayesian calibration framework and follows the statistical approach described in [Higdon et. al. (2004)](http://epubs.siam.org/doi/abs/10.1137/S1064827503426693).

Stan is portable across [many computing environments](http://mc-stan.org/users/interfaces/). Two interfaces, R ([main.R](https://github.com/adChong/bc-stan/blob/master/src/main.R)) and Python ([main.py](https://github.com/adChong/bc-stan/blob/master/src/main.py)) are provided here to interface with the stan models in this repository.

### Prerequisites

Install Stan and its required dependencies.
* RStan (R interface to stan): https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
* PyStan (Python interface to stan): http://pystan.readthedocs.io/en/latest/installation_beginner.html

### Usage

Main files
1. [main.R](https://github.com/adChong/bc-stan/blob/master/src/main.R): R interface for running Stan model
2. [main.py](https://github.com/adChong/bc-stan/blob/master/src/main.py): Python interface for running Stan model
3. [bcWithPred.stan](https://github.com/adChong/bc-stan/blob/master/src/bcWithPred.stan): Stan model for Bayesian calibration of building energy models with predictive inferences.
4. [bcWithoutPred.stan](https://github.com/adChong/bc-stan/blob/master/src/bcWithoutPred.stan): Stan model for Bayesian calibration of building energy models without predictive inferences.

To run Bayesian calibration with predictive inference in Stan, run [main.R](https://github.com/adChong/bc-stan/blob/master/src/main.R) or [main.py](https://github.com/adChong/bc-stan/blob/master/src/main.py) as is.

To run Bayesian calibration with predictive inference outside of Stan, comment lines 60-63 and line 85 and uncomment lines 65-68 and line 86 in [main.R](https://github.com/adChong/bc-stan/blob/master/src/main.R).

## Parameter screening using Morris method

The Python code for parameter screening is based on [Morris (1991)](http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00294) and implemented using the SALib python library.

### Usage

Main files
1. [sensitivity.py](https://github.com/adChong/bc-stan/blob/master/src/sensitivity.py): Python class for Morris method with E+ idf
2. [idf_functions.py](https://github.com/adChong/bc-stan/blob/master/src/idf_functions.py): Python class containing functions to modify E+ idf


### Prerequisites

Install eppy and SALib and their respective dependencies.
* eppy (scripting language for E+ idf files): https://pythonhosted.org/eppy/index.html
* SALib (Python library containing commonly used sensitivity analysis methods): https://salib.readthedocs.io

## Contact

If you need to get in touch for information about the code or its usage, you may drop us an [email](mailto:bdgczma@nus.edu.sg).






