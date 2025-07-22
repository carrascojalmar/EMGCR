compute_theta <- function(w, eta, link, tau) {
  if (link == "logit") {
    return(1 / (1 + exp(-w %*% eta)))
  } else if (link == "probit") {
    return(pnorm(w %*% eta))
  } else if (link == "plogit") {
    return((1 / (1 + exp(-w %*% eta)))^tau)
  } else if (link == "rplogit") {
    return(1 - (1 / (1 + exp(w %*% eta)))^tau)
  } else if (link == "cauchit") {
    return((1 / pi) * atan(w %*% eta) + 0.5)
  } else {
    stop("Unsupported link function")
  }
}
#' Plot multiple MCR model fits against Kaplan-Meier curve
#'
#' @import Formula
#' @import survival
#' @import knitr
#' @import flexsurv
#' @import tibble
#' @import stats
#' @importFrom actuar dinvgauss pinvgauss
#' @importFrom ggplot2 ggplot geom_step geom_line aes labs theme_minimal scale_y_continuous
#'
#' @param ... One or more fitted MCR objects from \code{MCRfit()}.
#'
#' @return A ggplot object with Kaplan-Meier and survival curves for each model.
#' @export


plot.MCR <- function(...) {
  fits <- list(...)

  if (length(fits) == 0) stop("At least one MCRfit object must be provided.")

  if (!all(sapply(fits, function(obj) inherits(obj, "MCR")))) {
    stop("All inputs must be objects of class 'MCR' from MCRfit().")
  }

  data <- fits[[1]]$data
  formula <- fits[[1]]$formula
  mf <- model.frame(Formula(formula), data = data)
  model.aux <- model.response(mf)
  y <- model.aux[, "time"]
  status <- model.aux[, "status"]
  times <- sort(unique(y))

  df_list <- list()

  for (i in seq_along(fits)) {
    fit <- fits[[i]]
    dist <- fit$dist
    formula <- fit$formula
    alpha <- fit$scale
    beta <- fit$coefficients
    eta <- fit$coefficients_cure
    link <- fit$link
    tau <- fit$tau
    data <- fit$data

    x <- model.matrix(Formula(formula), data = data, rhs = 1)
    w <- model.matrix(Formula(formula), data = data, rhs = 2)

    lambda <- exp(x %*% beta)
    theta <- compute_theta(w, eta, link, tau)

    sFit <- numeric(length(times))
    for (j in seq_along(times)) {
      t <- times[j]

      if (dist == "exponential") {
        surv_i <- pexp(t, rate = lambda, lower.tail = FALSE)
      } else if (dist == "rayleigh") {
        scale <- lambda^(-1 / 2)
        surv_i <- pweibull(t, shape = 2, scale = scale, lower.tail = FALSE)
      } else if (dist == "weibull") {
        scale <- lambda^(-1 / alpha)
        surv_i <- pweibull(t, shape = alpha, scale = scale, lower.tail = FALSE)
      } else if (dist == "lognormal") {
        surv_i <- plnorm(t, meanlog = -log(lambda), sdlog = alpha, lower.tail = FALSE)
      } else if (dist == "loglogistic") {
        surv_i <- flexsurv::pllogis(t, shape = alpha, scale = lambda, lower.tail = FALSE)
      } else if (dist == "invgauss") {
        surv_i <- actuar::pinvgauss(t, mean = lambda, shape = alpha, lower.tail = FALSE)
      } else {
        stop("Unsupported distribution: ", dist)
      }

      sFit[j] <- mean((1 - theta) + theta * surv_i)
    }

    df_list[[i]] <- data.frame(time = times, sFit = sFit, dist = dist)
  }

  df_all <- do.call(rbind, df_list)

  km_fit <- survfit(Surv(y, status) ~ 1, data = fits[[1]]$data)

  ggplot() +
    geom_step(aes(x = km_fit$time, y = km_fit$surv), color = "black",
              linetype = "dashed", size = 1, alpha = 0.8) +
    geom_line(data = df_all, aes(x = time, y = sFit,
                                 color = dist), size = 2)+
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Kaplan-Meier vs Multiple Fitted Distributions",
      x = "Time",
      y = "Survival Probability",
      color = "Distribution"
    ) +
    theme_minimal()
}
