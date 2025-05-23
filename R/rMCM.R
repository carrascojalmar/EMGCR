#' Generate Random Samples for Mixture Rate Cure (MRC) Model
#'
#' This function generates random samples for a specified MRC model, incorporating censoring and various distributional parameters.
#'
#' @param n An integer. The number of samples to generate.
#' @param x A numeric vector or matrix. The covariate(s) for the model.
#' @param w A numeric vector. The weights for each sample.
#' @param censor A numeric vector. Censoring values (e.g., limits for truncation or censoring).
#' @param beta A numeric vector. Coefficients for the covariates.
#' @param eta A numeric value. A scale parameter for the model.
#' @param alpha A numeric value. A shape parameter for the model.
#' @param link A character string. The link function to use. Defaults to `"logit"`. Options include `"logit"`, `"probit"`, `"rplogit"` and `"cauchit"`.
#' @param dist A character string. The distribution function to use. Defaults to `"weibull"`. Options include `"lognormal"` and `"loglogistic"`.
#' @param tau A numeric value. Threshold or cutoff parameter. Defaults to `1`.
#'
#' @return A data frame with the generated random samples, including:
#'   \item{sample}{The generated sample values.}
#'   \item{censoring}{The applied censoring values.}
#'
#' @details
#' This function simulates data based on a user-defined model. The link function transforms the covariates, while the parameters \code{beta}, \code{eta}, and \code{alpha} define the underlying distribution.
#'
#' @examples
#' # Example 1: Basic usage with default link function
#' n <- 500
#' alpha<-1.5
#' beta <- c(1,1,-2)
#' eta <- c(0.5,-0.6)
#'
#' p <- length(beta)
#' q <- length(eta)
#'
#' x <- matrix(rnorm(n*(p-1)),n,p-1)
#' x <- cbind(1,x)
#'
#' w <- matrix(runif(n*(q),-1,1),n,q) # no intercept for logit!
#' w <- scale(w)
#'
#' censoring_time<-10
#'
#' # sample
#' dataS<-rMCM(n, x, w, censor=censoring_time ,
#'            beta=beta, eta=eta, alpha=alpha)
#' dataAna0 <-data.frame(dataS$time,dataS$status,dataS$x[,-1],dataS$w)
#' colnames(dataAna0) <- c("time","status","x1","x2","w1","w2")
#' head(dataAna0,n=10)
#' @export
rMCM <- function(n,x,w,censor,alpha,beta,eta,
                 dist="weibull",
                 link="logit",tau=1){

  censoring_time <- censor
  # Probability of being *uncured* (eventually susceptible)

  if(link=="logit"){
    p_uncure <- 1 / (1 + exp(-w%*%eta))
  }else if(link=="probit"){
      p_uncure <- pnorm(w%*%eta)
  }else if(link=="plogit"){
        p_uncure <-(1/(1+exp(-w %*% eta)))^tau
  }else if(link=="rplogit"){
          p_uncure <-1-(1/(1+exp(w %*% eta)))^tau
  } else if(link=="cauchit"){
            p_uncure <-(1/pi)*atan(w%*%eta)+0.5
  }else{
    stop("Unsupported link function.")
  }

  uncured <- rbinom(n, 1, p_uncure)
  lambda <- exp(x %*% beta)
  event_time <- rep(NA, n)
  non_cured <- which(uncured == 1)

  u <- runif(length(non_cured))

  if(dist=="exponential"){
    event_time[non_cured] <- (-log(1 - u) / lambda[non_cured])^(1/1)
  }else if(dist=="rayleigh"){
    event_time[non_cured] <- (-log(1 - u)/lambda[non_cured])^(1/2)
  }else if(dist=="weibull"){
    event_time[non_cured] <- (-log(1 - u) / lambda[non_cured])^(1/alpha)
  }else if(dist=="lognormal"){
      event_time[non_cured] <- exp(qnorm(u)*alpha-log(lambda[non_cured]))
  }else if(dist == "loglogistic"){
        event_time[non_cured] <- lambda[non_cured]*(u/(1-u))^(1/alpha)
        }

  censor_time <- runif(n, 0, censoring_time)

  observed_time <- pmin(event_time, censor_time, na.rm = TRUE)

  # 0 = censored, 1 = event
  status <- ifelse(uncured == 0, 0, ifelse(event_time <= censor_time, 1, 0))


  #

  prop.cure.censur <-100*(sum(uncured==0)/n)
  df.Aux1 <- data.frame(uncured=uncured,status=status)
  df.Aux2 <- dplyr::filter(df.Aux1,uncured==1)
  prop.uncure.censur <- 100*(sum(df.Aux2$status==0)/dim(df.Aux2)[1])

  return(list(time=observed_time,status=status,x=x,w= w,
              pCcensur=prop.cure.censur,
              pUCcensur=prop.uncure.censur))
}
