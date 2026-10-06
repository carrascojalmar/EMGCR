# EMGCR 0.3.0

## Breaking changes

* The latency component is now parametrized through the scale parameter
  `lambda = exp(x %*% beta)` for all distributions, as in an accelerated failure
  time model (the same convention as `survival::survreg()`). Previously, the
  exponential, Rayleigh, Weibull and log-normal distributions used a rate
  parametrization. The fitted model and the log-likelihood are unchanged, but
  the latency coefficients of these four distributions now have the opposite
  sign (and, for Weibull and Rayleigh, a different scale) relative to version
  0.2.0. A positive coefficient now means longer survival for the uncured. The
  log-logistic and inverse Gaussian results are unchanged. `rMCM()`,
  `residuals()`, `qqMCR()` and `plot()` follow the new parametrization.
* The `liver` dataset was replaced by the liver cancer data used in the
  accompanying paper (1,736 patients; covariates `sex`, `age`, `medh`, `race`,
  `grade`, `chemo` and `radio`).

## Bug fixes

* Standard errors of the latency coefficients for the log-normal distribution
  were too small; they are now correct.
* `summary()` failed when the data had been modified after the fit (for example,
  when running the help-page examples).
* `rMCM()` failed when `x` had a single covariate besides the intercept.
* `plot()` no longer triggers the ggplot2 deprecation warning about `size`.
* For the inverse Gaussian distribution, the starting values of the latency
  coefficients followed the old rate parametrization, which could make the EM
  algorithm stop at a degenerate solution for some models.

## New features and improvements

* `summary()` shows significance stars and a summary of the quantile residuals,
  and returns the coefficient tables as numeric matrices.
* `residuals()` returns an object that prints the first residuals and a summary
  instead of every value.
* `MCRfit()` warns when the EM algorithm does not converge within `maxit`
  iterations and stores a logical `convergence` element; `summary()` reports it.
* `print()` and `summary()` show the link function (with `tau` only for the
  power links), label the incidence coefficients consistently as the
  "Uncured part", show the shape parameter only when it is estimated, and
  return their argument invisibly.
* Documentation of `MCRfit()`, `rMCM()` and `plot()` was corrected and expanded.
* `knitr` is no longer a dependency.
