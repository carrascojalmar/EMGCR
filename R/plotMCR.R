#' Compare fitted MCR models with the Kaplan-Meier curve
#'
#' Plots the Kaplan-Meier estimate of the survival function together with the
#' population survival function implied by one or more fitted mixture cure rate
#' models.
#'
#' @details
#' For each model, the fitted population survival function is averaged over the
#' observations,
#' \deqn{\hat S_{pop}(t) = \frac{1}{n} \sum_{i=1}^n \{1 - \hat\theta_i + \hat\theta_i \hat S(t \mid x_i)\},}
#' where \eqn{\hat\theta_i} is the estimated probability of being uncured and
#' \eqn{\hat S(t \mid x_i)} is the estimated survival function of the uncured.
#' The Kaplan-Meier curve (dashed) is computed from the data of the first model,
#' so all models should be fitted to the same data. Curves are coloured by
#' latency distribution.
#'
#' @import Formula
#' @import survival
#' @import flexsurv
#' @import tibble
#' @import stats
#' @importFrom actuar dinvgauss pinvgauss
#' @importFrom ggplot2 ggplot geom_step geom_line aes labs theme_minimal scale_y_continuous
#'
#' @param ... One or more fitted MCR objects from \code{MCRfit()}.
#'
#' @return A \code{ggplot} object.
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
        surv_i <- pexp(t, rate = 1 / lambda, lower.tail = FALSE)
      } else if (dist == "rayleigh") {
        surv_i <- pweibull(t, shape = 2, scale = lambda, lower.tail = FALSE)
      } else if (dist == "weibull") {
        surv_i <- pweibull(t, shape = alpha, scale = lambda, lower.tail = FALSE)
      } else if (dist == "lognormal") {
        surv_i <- plnorm(t, meanlog = log(lambda), sdlog = alpha, lower.tail = FALSE)
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
              linetype = "dashed", linewidth = 1, alpha = 0.8) +
    geom_line(data = df_all, aes(x = time, y = sFit,
                                 color = dist), linewidth = 2)+
    scale_y_continuous(limits = c(0, 1)) +
    labs(
      title = "Kaplan-Meier vs Multiple Fitted Distributions",
      x = "Time",
      y = "Survival Probability",
      color = "Distribution"
    ) +
    theme_minimal()
}
