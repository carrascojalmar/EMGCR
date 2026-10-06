#' Compute residuals for MCR model
#'
#' This function computes Global Cox-Snell and randomized quantile residuals for objects of class \code{MCR}.
#'
#' @import Formula
#' @importFrom survival Surv
#' @import flexsurv
#' @import tibble
#' @import stats
#' @importFrom actuar dinvgauss pinvgauss
#'
#' @param object An object of class \code{MCR}, typically returned from \code{\link{MCRfit}}.
#' @param type Type of residual.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of residuals of class \code{"MCRresiduals"}. Printing it shows
#'   the first six residuals and a summary; use \code{as.vector()} to drop the class.
#' @export
#' @method residuals MCR
#'
#' @examples
#' data(liver)
#' liver$sex <- factor(liver$sex)
#' liver$grade <- factor(liver$grade)
#' liver$radio <- factor(liver$radio)
#' liver$chemo <- factor(liver$chemo)
#'
#' model <- MCRfit(
#'   survival::Surv(time, status) ~ age + sex + grade + radio + chemo |
#'     age + medh + grade + radio + chemo,
#'   dist = "loglogistic", link = "plogit", tau = 0.15,
#'   data = liver
#' )
#' residuals(model, type = "quantile")
#'
residuals.MCR <- function(object, type = c("cox-snell","quantile"), ...) {

  type <- match.arg(type)

  formula <- object$formula
  data <-object$data
  n <- object$n

  mf1 <- model.frame(Formula(formula), data = data)
  x <- model.matrix(Formula(formula), data = mf1, rhs = 1)
  w <- model.matrix(Formula(formula), data = mf1, rhs = 2)

  mf <- model.frame(Formula(formula), data = data)
  model.aux <- model.response(mf)
  cc <- model.aux[, "status"]
  y <- model.aux[, "time"]

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
    aux <- pexp(q=y,rate=1/lambda,lower.tail = FALSE)
  } else if (dist == "rayleigh") {
    aux <- pweibull(q=y,shape=2,scale=lambda,lower.tail = FALSE)
  } else if (dist == "weibull") {
    aux <- pweibull(q=y,shape=alpha,scale=lambda,
             lower.tail = FALSE)
  } else if (dist == "lognormal") {
    aux <- plnorm(q=y,meanlog = log(lambda),
           sdlog = alpha,lower.tail = FALSE)
  } else if (dist == "loglogistic") {
    aux <- flexsurv::pllogis(q=y, shape=alpha,scale=lambda,
                      lower.tail = FALSE)
  } else if (dist == "invgauss"){
    aux <- actuar::pinvgauss(y, mean = lambda, shape = alpha, lower.tail = FALSE)
  } else {
    stop("Unsupported distribution")
  }

  auxR <- theta*aux+(1-theta)

  if (type == "cox-snell") {
    resCM <- -log(auxR)
  } else if (type == "quantile") {
    resCM <- qnorm(cc*(1-auxR)+(1-cc)*runif(n,1-auxR))
  }

  resCM <- as.vector(resCM)
  names(resCM) <- rownames(mf)
  structure(resCM, type = type, class = "MCRresiduals")
}

#' @export
print.MCRresiduals <- function(x, n = 6L, digits = max(3L, getOption("digits") - 3L), ...) {
  label <- if (identical(attr(x, "type"), "quantile")) "Quantile" else "Cox-Snell"
  res <- as.vector(x)
  names(res) <- names(x)

  cat(label, " residuals (n = ", length(res), "):\n", sep = "")
  print(res[seq_len(min(n, length(res)))], digits = digits)
  if (length(res) > n) cat("... (", length(res) - n, " more)\n", sep = "")

  cat("\nSummary:\n")
  print(summary(res), digits = digits)
  invisible(x)
}
