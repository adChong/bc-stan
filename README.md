# bc-stan

bc-stan provides code for the Bayesian calibration of Building Energy Models using the Stan modeling language. The stan code is based on [Kennedy and O'Hagan's (2001)](http://onlinelibrary.wiley.com/doi/10.1111/1467-9868.00294/abstract) Bayesian calibration framework and follows the statistical approach described in [Higdon et. al. (2004)](http://epubs.siam.org/doi/abs/10.1137/S1064827503426693).

Stan is portable across [many computing environments](http://mc-stan.org/users/interfaces/). Two interfaces, R ([main.R](https://github.com/adChong/bc-stan/blob/master/src/main.R)) and Python (main.py) are provided here to interface with the stan models in this repository.

## Table of Contents

1. [Prerequisits](#Prerequisites)
2. [Usage](#Usage)
3. [Documentation](#Documentation)
4. [Contact](#Contact)

## Prerequisites

Install Stan and its required dependencies.
* RStan (R interface to stan): https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
* PyStan (Python interface to stan): http://pystan.readthedocs.io/en/latest/installation_beginner.html

## Usage

Main files
1. [main.R](https://github.com/adChong/bc-stan/blob/master/src/main.R): R interface for running Stan model
2. main.py: Python interface for running Stan model
3. [bcWithPred.stan](https://github.com/adChong/bc-stan/blob/master/src/bcWithPred.R): Stan model for Bayesian calibration of building energy models with predictive inferences.
4. [bcWithoutPred.stan](https://github.com/adChong/bc-stan/blob/master/src/bcWithoutPred.R): Stan model for Bayesian calibration of building energy models without predictive inferences.



## Documentation

## Contact

If you need to get in touch for information about the code or its usage, you may drop us an [email](bdgczma@nus.edu.sg).






