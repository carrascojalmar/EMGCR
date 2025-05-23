#' Compute residuals for MCR model
#'
#' This function computes Global Cox-Snell and randomized quantile residuals for objects of class \code{MCR}.
#'
#' @param object An object of class \code{MCR}, typically returned from \code{\link{MCRfit}}.
#' @param type Type of residual.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of residuals.
#' @export
#' @method residuals MCR
#'
#' @examples
#' data(liver)
#' #head(liver)
#' model <- MCRfit(Surv(time,status)~age+medh+relapse+grade|-1+sex+age+medh+grade,data=liver)
#' summary(residuals(model,type="quantile"))
#'
residuals.MCR <- function(object, type = c("cox-snell","quantile"), ...) {

  type <- match.arg(type)

  formula <- object$formula
  data <-object$data
  n <- object$n

  mf <- model.frame(Formula(formula), data = data)
  model.aux <- model.response(mf)
  cc <- model.aux[, "status"]
  y <- model.aux[, "time"]

  x <- model.matrix(Formula(formula), data = data, rhs = 1)
  w <- model.matrix(Formula(formula), data = data, rhs = 2)

  beta <- object$coefficients
  eta <- object$coefficients_cure
  alpha <- object$scale
  dist <- object$dist
  link <- object$link
  tau <- object$tau

  lambda <- exp(x %*% beta)

  if (link == "logit") {
    theta <- 1 / (1 + exp(-w %*% eta))
  } else if (link == "probit") {
    theta <- pnorm(w %*% eta)
  } else if (link == "plogit") {
    theta <- (1 / (1 + exp(-w %*% eta)))^tau
  } else if (link == "rplogit") {
    theta <- 1 - (1 / (1 + exp(w %*% eta)))^tau
  } else if (link == "cauchit") {
    theta <- (1 / pi) * atan(w %*% eta) + 0.5
  } else {
    stop("Unsupported link function")
  }

  if (dist == "exponential") {
    aux <- exponential_sf(y, alpha = 1, lambda)
  } else if (dist == "rayleigh") {
    aux <- rayleigh_sf(y, alpha = 2, lambda)
  } else if (dist == "weibull") {
    aux <- weibull_sf(y, alpha, lambda)
  } else if (dist == "lognormal") {
    aux <- lognormal_sf(y, alpha, -log(lambda))
  } else if (dist == "loglogistic") {
    aux <- loglogistic_sf(y, alpha, lambda)
  } else {
    stop("Unsupported distribution")
  }

  auxR <- theta*aux+(1-theta)

  if (type == "cox-snell") {
    resCM <- -log(auxR)
    return(resCM)

  } else if (type == "quantile") {
    resCM <- qnorm(cc*(1-auxR)+(1-cc)*runif(n,1-auxR))
    return(resCM)
  }

}
