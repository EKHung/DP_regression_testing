# Minimax signal detection in regression under differential privacy

This repository contains the R code for the numerical experiments in the paper. 

## Setup

You will need R, the packages `Rcpp`, `ggplot2` (version 3.4.0 or later) and `xtable`, and a working C++ compiler. 
Run the scripts from the directory containing the files below. 

## Files

| File | Description |
| --- | --- |
| `test_functions.R` | Implements the transformed-mean, coordinate and response-only tests, together with the comparison methods and shared helper functions. |
| `run_gaussian.R` | Compares five private tests under Gaussian linear regression. |
| `run_illustrations.R` | Runs the Cauchy and exponential examples, the response-only dimension experiment, and the sparsity experiment. |
| `run_sensitivity.R` | Examines sensitivity to tuning constants, kernel bandwidths and aggregation block counts. |
| `plot_regression_results.R` | Produces the figures and LaTeX tables from the saved simulation results. |


## Running the experiments

From a terminal in this directory, run:

```sh
Rscript run_gaussian.R
Rscript run_illustrations.R
Rscript run_sensitivity.R
```
These scripts save `gauss.RData`, `illustrations.RData` and `sensitivity.Rdata`, respectively. 

Once the three result files are available, run:

```sh
Rscript plot_regression_results.R
```

The plotting script creates five figures in the working directory:

- `comparison.png`: the Gaussian comparison;
- `cauchy.png`: the Cauchy example;
- `exponential.png`: the exponential example;
- `dimension_stability.png`: the response-only test across dimensions;
- `coordinate_sparsity.png`: the coordinate test across sparsity levels.

Error bars show one Monte Carlo standard error on either side of the estimated sum of type-I and type-II errors. The script also prints four LaTeX sensitivity tables to the console. The tables use the `booktabs` package.

