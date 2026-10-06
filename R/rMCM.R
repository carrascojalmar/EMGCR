#' Generate Random Samples for Mixture Cure Rate (MCR) Model
#'
#' Simulates survival data from a mixture cure rate model with covariates, a chosen link function for the incidence part and a chosen latency distribution. Censoring times are drawn from a uniform distribution on \eqn{(0, \code{censor})}.
#'
#' @import Formula
#' @importFrom survival Surv
#' @import flexsurv
#' @import tibble
#' @import stats
#' @importFrom actuar dinvgauss pinvgauss rinvgauss
#'
#' @param n Integer. Number of observations to simulate.
#' @param x Matrix or numeric. Covariate matrix for the latency component. The first column must be the intercept (a column of ones).
#' @param w Matrix or numeric. Covariate matrix for the incidence part (the probability of being uncured). Include a column of ones if an intercept is wanted.
#' @param censor Numeric. Maximum censoring time (uniformly distributed).
#' @param alpha Numeric. Shape parameter \eqn{\alpha} of the latency distribution. Ignored for the exponential and Rayleigh distributions.
#' @param beta Numeric vector. Coefficients for the latency part, which enter through the scale parameter \eqn{\lambda = \exp(x^\top \beta)} (see \code{\link{MCRfit}} for details).
#' @param eta Numeric vector. Coefficients for the incidence part.
#' @param dist Character. Distribution for the latency part. Options: `"weibull"`, `"lognormal"`, `"loglogistic"`, `"invgauss"`, `"exponential"`, `"rayleigh"`.
#' @param link Character. Link function for the probability of being uncured. Options: `"logit"`, `"probit"`,`"plogit"` ,`"rplogit"`, `"cauchit"`.
#' @param tau A numeric value used when \code{link = "plogit"} or \code{"rplogit"}. Defaults to 1.
#'
#' @return A \code{\link[tibble]{tibble}} with columns:
#' \describe{
#'   \item{time}{Observed (possibly censored) survival time.}
#'   \item{status}{Event indicator (1 = event, 0 = censored).}
#'   \item{x1, x2, ...}{Columns of \code{x} without the intercept.}
#'   \item{w1, w2, ...}{Columns of \code{w}.}
#' }
#' It also has two attributes: \code{"pCcensur"}, the percentage of cured
#' individuals, and \code{"pUCcensur"}, the percentage of censored observations
#' among the uncured.
#'
#' @examples
#' # Example: Simulating survival data using the inverse Gaussian distribution
#' library(EMGCR)
#'
#' n <- 500
#' beta <- c(1, -1, -2)
#' eta <- c(0.5, -0.5)
#' alpha <- 1.5
#'
#' p <- length(beta)
#' q <- length(eta)
#'
#' set.seed(10)
#' X <- matrix(rnorm(n*(p-1),0,1),n,p-1)
#' X <- cbind(1,X)
#'
#' set.seed(20)
#' W <- matrix(runif(n*q,-1,1),n,q)
#' W <- scale(W)
#'
#' max_censoring <- 10
#'
#' set.seed(1234)
#' sim_data <- rMCM(n=n, x = X, w = W,
#'                  censor = max_censoring,
#'                  beta = beta, eta = eta,
#'                  alpha = alpha,
#'                  link = "logit", dist = "invgauss", tau = 1)
#'
#' names(sim_data)
#' head(sim_data)
#' attributes(sim_data)
#' attr(sim_data, "pCcensur")
#' attr(sim_data, "pUCcensur")
#' @export
rMCM <- function(n, x, w, censor, alpha, beta, eta,
                 dist = "weibull", link = "logit", tau = 1) {

  # Input checks
  stopifnot(is.numeric(n), n > 0,
            is.matrix(x) || is.numeric(x),
            is.matrix(w) || is.numeric(w),
            length(beta) == ncol(x),
            length(eta) == ncol(w),
            dist %in% c("weibull", "lognormal", "loglogistic", "invgauss", "exponential", "rayleigh"),
            link %in% c("logit", "probit", "plogit", "rplogit", "cauchit"))

  x <- as.matrix(x)
  w <- as.matrix(w)
  lambda <- exp(x %*% beta)

  # Compute uncure probabilities
  linpred <- w %*% eta
  p_uncure <- switch(link,
                     logit = 1 / (1 + exp(-linpred)),
                     probit = pnorm(linpred),
                     plogit = (1 / (1 + exp(-linpred)))^tau,
                     rplogit = 1 - (1 / (1 + exp(linpred)))^tau,
                     cauchit = (1 / pi) * atan(linpred) + 0.5)

  # Cure status and initialize times
  uncured <- rbinom(n, 1, p_uncure)
  event_time <- rep(NA_real_, n)
  idx <- which(uncured == 1)
  u <- runif(length(idx))

  # Latency model
  if (dist == "exponential") {
    event_time[idx] <- (-log(1 - u)) * lambda[idx]
  } else if (dist == "rayleigh") {
    event_time[idx] <- (-log(1 - u))^(1 / 2) * lambda[idx]
  } else if (dist == "weibull") {
    event_time[idx] <- (-log(1 - u))^(1 / alpha) * lambda[idx]
  } else if (dist == "lognormal") {
    event_time[idx] <- exp(qnorm(u) * alpha + log(lambda[idx]))
  } else if (dist == "loglogistic") {
    event_time[idx] <- lambda[idx] * (u / (1 - u))^(1 / alpha)
  } else if (dist == "invgauss") {
    event_time[idx] <- rinvgauss(length(idx), lambda[idx], shape = alpha)
  }

  censor_time <- runif(n, 0, censor)
  observed_time <- pmin(event_time, censor_time, na.rm = TRUE)
  status <- ifelse(uncured == 0, 0, as.numeric(event_time <= censor_time))

  # Proportions
  pCcensur  <- 100 * mean(uncured == 0)
  pUCcensur <- ifelse(sum(uncured == 1) > 0,
                      100 * mean(status[uncured == 1] == 0),
                      NA_real_)

  xx <- x[, -1, drop = FALSE]
  colnames(xx) <- paste0("x", seq_len(ncol(x)-1))
  colnames(w) <- paste0("w", seq_len(ncol(w)))

  out <- tibble::tibble(
    time = observed_time,
    status = status,
    !!!as.data.frame(xx),
    !!!as.data.frame(w)
  )

  attr(out, "pCcensur")  <- pCcensur
  attr(out, "pUCcensur") <- pUCcensur

  return(out)
}
