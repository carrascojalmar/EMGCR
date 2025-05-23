envelopeRQR <- function(x, nsim = 100) {
  n <- length(x)
  x_sorted <- sort(x)
  x_theo <- qnorm(ppoints(n))

  sim_res <- matrix(rnorm(nsim * n), nrow = nsim, ncol = n)
  sim_sorted <- t(apply(sim_res, 1, sort))


  lower <- apply(sim_sorted, 2, quantile, probs = 0.025)
  upper <- apply(sim_sorted, 2, quantile, probs = 0.975)
  mean_env <- apply(sim_sorted, 2, mean)


  y_range <- range(x_sorted, lower, upper)

  # Gráfico
  plot(x_theo, x_sorted,
       xlab = "Theoretical Quantiles (N(0,1))",
       ylab = "Quantile Residuals",
       main = "Normal QQ-Plot of Quantile Residuals with Envelope",
       pch = 19, col = "steelblue",
       ylim = y_range)
  lines(x_theo, lower, col = "gray", lty = 2, lwd = 1.2)
  lines(x_theo, mean_env, col = "black", lwd = 2)
  lines(x_theo, upper, col = "gray", lty = 2, lwd = 1.2)
}
#' QQ-Plot of Residuals for MCR Model
#'
#' @param object An object of class "MCR"
#' @param type Type of residual
#' #' @param envelope Logical; if TRUE, the envelope simulation will be performed.
#' @param ... Additional arguments (not used)
#' @examples
#' data(liver)
#' #head(liver)
#' model <- MCRfit(Surv(time,status)~age+medh+relapse+grade|-1+sex+age+medh+grade,data=liver)
#' qqMCR(model,type="quantile",envelope=TRUE)
#
#' @export
qqMCR <- function(object, type = c("cox-snell", "quantile"),
                  envelope = FALSE, nsim = 100, censor = NULL, ...) {

  type <- match.arg(type)

  data    <- object$data
  n       <- object$n
  formula <- object$formula
  x       <- model.matrix(Formula(formula), data = data, rhs = 1)
  w       <- model.matrix(Formula(formula), data = data, rhs = 2)

  beta  <- object$coefficients
  eta   <- object$coefficients_cure
  alpha <- object$scale
  dist  <- object$dist
  link  <- object$link
  tau   <- object$tau

  res <- residuals(object, type)
  res <- res[is.finite(res)]

  if (!envelope) {
    switch(type,
           "cox-snell" = {
             theo_q <- qexp(ppoints(length(res)))
             qqplot(theo_q, sort(res),
                    main = "Exponential QQ-Plot of Residuals",
                    xlab = "Theoretical Quantiles (Exp(1))",
                    ylab = "Empirical Residuals",
                    pch = 19, col = "steelblue")
           },
           "quantile" = {
             qqnorm(res,
                    main = "Normal QQ-Plot of Residuals",
                    xlab = "Theoretical Quantiles (N(0,1))",
                    ylab = "Quantile Residuals",
                    pch = 19, col = "steelblue")
           })
    return(invisible())
  }

  if (type == "cox-snell") {

    if (is.null(censor)) stop("Argument 'censor' must be provided when envelope = TRUE and type = 'cox-snell'.")

    mres <- matrix(NA, nrow = n, ncol = nsim)

    for (r in 1:nsim) {
      df.aux <- rMCM(n, x, w, censor = censor, beta = beta, eta = eta,
                     alpha = alpha, link = link, dist = dist, tau = tau)

      data$time   <- df.aux$time
      data$status <- df.aux$status

      model.env <- suppressWarnings(
        MCRfit(formula, link = link, dist = dist, tau = tau, data = data)
      )

      res.env <- residuals(model.env, type)
      mres[, r] <- sort(res.env)
    }

    min   <- apply(mres, 1, min)
    med   <- apply(mres, 1, mean)
    max   <- apply(mres, 1, max)
    theo_q <- qexp(ppoints(n))

    plot(theo_q, sort(res),
         main = "Exponential QQ-Plot of Residuals with Envelope",
         xlab = "Theoretical Quantiles (Exp(1))",
         ylab = "Empirical Residuals",
         pch = 19, col = "steelblue")
    lines(theo_q, min, col = "gray", lty = 2, lwd = 2)
    lines(theo_q, med, col = "black", lwd = 2)
    lines(theo_q, max, col = "gray", lty = 2, lwd = 2)

  } else if (type == "quantile") {
    envelopeRQR(res, nsim)
  }
}
