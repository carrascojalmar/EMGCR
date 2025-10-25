# link function

plogit <- function(tau) {

  # Link function: maps from probabilities (mu) to the linear predictor (eta)
  linkfun <- function(mu) {
    (log(mu^(1/tau) / (1 - mu^(1/tau)))) # Standard logit link
  }

  # Inverse link function: maps from the linear predictor (eta) to probabilities (mu)
  linkinv <- function(eta) {
    (1 / (1 + exp(-eta)))^tau  # Logistic function raised to power tau
  }

  # Derivative of the inverse link function with respect to eta
  mu.eta <- function(eta) {
    logistic <- 1 / (1 + exp(-eta))  # Standard logistic function
    tau * logistic^(tau - 1) * exp(-eta) / (1 + exp(-eta))^2  # Derivative with respect to eta
  }

  # Valid eta check (always valid for this transformation)
  valideta <- function(eta) TRUE

  # Create and return the custom link list that glm() or vglm() will use
  link <- list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
               valideta = valideta,
               name = paste("plogit(tau=", tau, ")",
                            sep=""))
  class(link) <- "link-glm"
  return(link)
}

#  Reversal Power Logit (type II)

rplogit <- function(tau) {

  # Link function: maps from probabilities (mu) to the linear predictor (eta)
  linkfun <- function(mu) {
    log((1-mu)^(-1/tau)-1) # Standard logit link
  }

  # Inverse link function: maps from the linear predictor (eta) to probabilities (mu)
  linkinv <- function(eta) {
    1-(1 / (1 + exp(eta)))^tau  # Logistic function raised to power tau
  }

  # Derivative of the inverse link function with respect to eta
  mu.eta <- function(eta) {
    logisticr <- 1 / (1 + exp(eta))  # Standard logistic function
    tau * logisticr^(tau+1) * exp(eta)  # Derivative with respect to eta
  }

  # Valid eta check (always valid for this transformation)
  valideta <- function(eta) TRUE

  # Create and return the custom link list that glm() or vglm() will use
  link <- list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, valideta = valideta, name = paste("rplogit(tau=", tau, ")", sep=""))
  class(link) <- "link-glm"
  return(link)
}

# likelihood
# Link function
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

# Log-likelihood calculation
likeMRC <- function(y, cc, x, w, alpha, beta, eta, link = "logit",
                    dist = "weibull", tau = 1) {
  lambda <- exp(x %*% beta)
  theta <- compute_theta(w, eta, link, tau)

  n <- length(y)
  ll <- numeric(n)

  for (i in 1:n) {
    time <- y[i]
    status <- cc[i]
    lam <- lambda[i]
    th <- theta[i]

    if (dist == "exponential") {
      dens <- dexp(time, rate = lam)
      surv <- pexp(time, rate = lam, lower.tail = FALSE)
    } else if (dist == "rayleigh") {
      scale <- lam^(-1/2)
      dens <- dweibull(time, shape = 2, scale = scale)
      surv <- pweibull(time, shape = 2, scale = scale, lower.tail = FALSE)
    } else if (dist == "weibull") {
      scale <- lam^(-1 / alpha)
      dens <- dweibull(time, shape = alpha, scale = scale)
      surv <- pweibull(time, shape = alpha, scale = scale, lower.tail = FALSE)
    } else if (dist == "lognormal") {
      dens <- dlnorm(time, meanlog = -log(lam), sdlog = alpha)
      surv <- plnorm(time, meanlog = -log(lam), sdlog = alpha, lower.tail = FALSE)
    } else if (dist == "loglogistic") {
      dens <- flexsurv::dllogis(time, shape = alpha, scale = lam)
      surv <- flexsurv::pllogis(time, shape = alpha, scale = lam, lower.tail = FALSE)
    } else if (dist == "invgauss") {
      dens <- actuar::dinvgauss(time, mean = lam, shape = alpha)
      surv <- actuar::pinvgauss(time, mean = lam, shape = alpha, lower.tail = FALSE)
    } else {
      stop("Unsupported distribution")
    }

    aux1 <- (th * dens)^status
    aux2 <- (1 - th + th * surv)^(1 - status)
    ll[i] <- aux1 * aux2
  }

  return(sum(log(pmax(ll, .Machine$double.xmin))))
}
#
gradInvGauss <- function(params, time, status,x,B){
  p <- ncol(x)
  beta <- as.vector(params[1:p])
  alpha <- params[p+1]
  lambda <- exp(x%*%beta)

  aux1 <- 1- pinvgauss(q=time, mean=lambda, shape=alpha,
                       lower.tail = TRUE, log.p = FALSE)

  if(length(which(aux1 == 0)) > 0) aux1[which(aux1 == 0)] <- .Machine$double.xmin

  aux2 <- (time/lambda)-1
  aux3 <- (time/lambda)+1

  z1 <- sqrt(alpha/time)*aux2
  z2 <- sqrt(alpha/time)*aux3


  U1 <- 0.5*((1/alpha)-((time-lambda)**2/(lambda*time)))
  U2 <- (aux2/2*sqrt(alpha*time))*dnorm(z1)+
    (2*exp(2*alpha/lambda)/lambda)*pnorm(-z2)-
    (aux3*exp(2*alpha/lambda)/sqrt(alpha*time))*dnorm(-z2)

  Ualpha <- sum(B*status*U1-(1-status)*(U2/aux1))

  Ubeta <- numeric()
  for (j in 1:p){
    xj <- x[,j]

    aux4 <- alpha*(time-lambda)*xj/time
    U3 <- aux4*(1+((time-lambda)/2*lambda))

    U4 <- -sqrt(alpha/time)*(time*xj/lambda)*dnorm(z1)-
      (2*alpha/lambda)*exp(2*alpha/lambda)*pnorm(-z2)+
      sqrt(alpha/time)*(time*xj/lambda)*exp(2*alpha/lambda)*dnorm(-z2)

    Ubeta[j] <- sum(B*status*U3-(1-status)*(U4/aux1))
  }

  U <- c(Ubeta,Ualpha)

  return(U)
}
likeInvGauss <- function(params, time, status, x, B){
  p <- ncol(x)
  beta <- as.vector(params[1:p])
  alpha <- params[p+1]
  lambda <- exp(x%*%beta)

  aux1<-B*status*dinvgauss(x=time, mean=lambda,
                           shape=alpha, log = TRUE)
  aux11<-pinvgauss(q=time, mean=lambda, shape=alpha,
                   lower.tail = FALSE, log.p = FALSE)

  if(length(which(aux11 == 0)) > 0) aux11[which(aux11 == 0)] <- .Machine$double.xmin
  if (any(is.nan(aux11) | is.infinite(aux11))) return(1e+10)

  aux2<-B*(1-status)*log(aux11)

  return(-sum((aux1+aux2)))
}
#
epMCR <- function(y, cc, x, w, B, alpha, beta, eta, tau,
                  dist = dist, link = link) {

  n <- length(y)
  p <- ncol(x)
  q <- ncol(w)

  if(dist == "exponential" && link=="probit"){
    alpha <- 1
    lambda<- exp(x %*% beta)
    theta<-pnorm(w %*% eta)
    pdfN<-dnorm(w %*% eta)
    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<- (B[i]*pdfN[i]/(theta[i])-(1-B[i])*pdfN[i]/(1-theta[i]))%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "exponential" && link=="logit"){

    alpha <- 1
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))

    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<- (B[i]*(1-theta[i])-(1-B[i])*theta[i])%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "exponential" && link=="plogit"){

    alpha <- 1
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))^tau
    derv<-tau*exp(-w %*% eta)/(1+exp(-w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "exponential" && link=="rplogit"){

    alpha <- 1
    lambda<- exp(x %*% beta)
    theta<-1-(1/(1+exp(w %*% eta)))^tau
    derv<-tau*exp(w %*% eta)/(1+exp(w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "exponential" && link=="cauchit"){

    alpha <- 1
    lambda<- exp(x %*% beta)
    theta<-(1/pi)*atan(w%*%eta)+0.5
    derv<-(1/pi)*(1/(1+(w%*%eta)^2))
    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "rayleigh" && link=="probit"){
    alpha <- 2
    lambda<- exp(x %*% beta)
    theta<-pnorm(w %*% eta)
    pdfN<-dnorm(w %*% eta)
    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<- (B[i]*pdfN[i]/(theta[i])-(1-B[i])*pdfN[i]/(1-theta[i]))%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "rayleigh" && link=="logit"){

    alpha <- 2
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))

    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<- (B[i]*(1-theta[i])-(1-B[i])*theta[i])%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "rayleigh" && link=="plogit"){

    alpha <- 2
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))^tau
    derv<-tau*exp(-w %*% eta)/(1+exp(-w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "rayleigh" && link=="rplogit"){

    alpha <- 2
    lambda<- exp(x %*% beta)
    theta<-1-(1/(1+exp(w %*% eta)))^tau
    derv<-tau*exp(w %*% eta)/(1+exp(w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "rayleigh" && link=="cauchit"){

    alpha <- 2
    lambda<- exp(x %*% beta)
    theta<-(1/pi)*atan(w%*%eta)+0.5
    derv<-(1/pi)*(1/(1+(w%*%eta)^2))
    MI<-matrix(0,(p+q),(p+q))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      #Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2),(p+q),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "weibull" && link=="probit"){
    lambda<- exp(x %*% beta)
    theta<-pnorm(w %*% eta)
    pdfN<-dnorm(w %*% eta)
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<- (B[i]*pdfN[i]/(theta[i])-(1-B[i])*pdfN[i]/(1-theta[i]))%*%(w[i,])
      Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2,Aux3),(p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "weibull" && link=="logit"){
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))

    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<- (B[i]*(1-theta[i])-(1-B[i])*theta[i])%*%(w[i,])
      Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2,Aux3),(p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "weibull" && link=="plogit"){
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))^tau
    derv<-tau*exp(-w %*% eta)/(1+exp(-w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2,Aux3),(p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "weibull" && link=="rplogit"){
    lambda<- exp(x %*% beta)
    theta<-1-(1/(1+exp(w %*% eta)))^tau
    derv<-tau*exp(w %*% eta)/(1+exp(w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2,Aux3),(p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }

  if(dist == "weibull" && link=="cauchit"){
    lambda<- exp(x %*% beta)
    theta<-(1/pi)*atan(w%*%eta)+0.5
    derv<-(1/pi)*(1/(1+(w%*%eta)^2))
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      Aux1<-(cc[i]*B[i]*(1-y[i]^alpha*lambda[i])-(1-cc[i])*B[i]*(y[i]^alpha*lambda[i]))%*%(x[i,])
      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      Aux3<-cc[i]*B[i]*(1/alpha+log(y[i])-lambda[i]*y[i]^alpha*log(y[i]))- (1-cc[i])*B[i]*(lambda[i]*y[i]^alpha*log(y[i]))
      Aux4<-matrix(c(Aux1,Aux2,Aux3),(p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)
    }
  }


  if(dist == "lognormal" && link=="probit"){
    mu <- (x %*% beta)
    theta<-pnorm(w %*% eta)
    pdfN<-dnorm(w %*% eta)

    Z <- (log(y) - mu)/alpha
    pdfZ <- dnorm(Z)
    cdfZ <- pnorm(Z)
    S_y <- 1 - cdfZ

    MI<-matrix(0,(p+q+1),(p+q+1))
    for (i in 1:n){
      S_y[i] <- max(S_y[i], .Machine$double.eps)
      Aux1 <- (cc[i]*B[i]*(Z[i]/alpha)+
                 (1-cc[i])*B[i]*(pdfZ[i]/(S_y[i]*alpha)))%*%(x[i,])
      Aux2<- (B[i]*pdfN[i]/(theta[i])-(1-B[i])*pdfN[i]/(1-theta[i]))%*%(w[i,])
      Aux3 <- cc[i]*B[i]*(-1/(alpha)+(Z[i]^2)/(alpha))+
        (1-cc[i])*B[i]*(Z[i]*pdfZ[i])/(alpha*S_y[i])
      Aux4 <- matrix(c(Aux1,Aux2,Aux3), (p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)

    }

  }

  if(dist == "lognormal" && link=="logit"){
    mu <- (x %*% beta)
    theta<-1/(1+exp(-w %*% eta))
    Z <- (log(y) - mu)/alpha
    pdfZ <- dnorm(Z)
    cdfZ <- pnorm(Z)
    S_y <- 1 - cdfZ

    MI<-matrix(0,(p+q+1),(p+q+1))
    for (i in 1:n){
      S_y[i] <- max(S_y[i], .Machine$double.eps)
      Aux1 <- (cc[i]*B[i]*(Z[i]/alpha)+
                 (1-cc[i])*B[i]*(pdfZ[i]/(S_y[i]*alpha)))%*%(x[i,])
      Aux2<- (B[i]*(1-theta[i])-(1-B[i])*theta[i])%*%(w[i,])
      Aux3 <- cc[i]*B[i]*(-1/(alpha)+(Z[i]^2)/(alpha))+
        (1-cc[i])*B[i]*(Z[i]*pdfZ[i])/(alpha*S_y[i])
      Aux4 <- matrix(c(Aux1,Aux2,Aux3), (p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)

    }

  }

  if(dist == "lognormal" && link=="plogit"){
    mu <- (x %*% beta)
    theta<-1/(1+exp(-w %*% eta))^tau
    derv<-tau*exp(-w %*% eta)/(1+exp(-w %*% eta))^(tau+1)

    Z <- (log(y) - mu)/alpha
    pdfZ <- dnorm(Z)
    cdfZ <- pnorm(Z)
    S_y <- 1 - cdfZ

    MI<-matrix(0,(p+q+1),(p+q+1))
    for (i in 1:n){
      S_y[i] <- max(S_y[i], .Machine$double.eps)
      Aux1 <- (cc[i]*B[i]*(Z[i]/alpha)+
                 (1-cc[i])*B[i]*(pdfZ[i]/(S_y[i]*alpha)))%*%(x[i,])
      Aux2<- (B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      Aux3 <- cc[i]*B[i]*(-1/(alpha)+(Z[i]^2)/(alpha))+
        (1-cc[i])*B[i]*(Z[i]*pdfZ[i])/(alpha*S_y[i])
      Aux4 <- matrix(c(Aux1,Aux2,Aux3), (p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)

    }

  }


  if(dist == "lognormal" && link=="rplogit"){
    mu <- (x %*% beta)
    theta<-1-(1/(1+exp(w %*% eta)))^tau
    derv<-tau*exp(w %*% eta)/(1+exp(w %*% eta))^(tau+1)

    Z <- (log(y) - mu)/alpha
    pdfZ <- dnorm(Z)
    cdfZ <- pnorm(Z)
    S_y <- 1 - cdfZ

    MI<-matrix(0,(p+q+1),(p+q+1))
    for (i in 1:n){
      S_y[i] <- max(S_y[i], .Machine$double.eps)
      Aux1 <- (cc[i]*B[i]*(Z[i]/alpha)+
                 (1-cc[i])*B[i]*(pdfZ[i]/(S_y[i]*alpha)))%*%(x[i,])
      Aux2<- (B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      Aux3 <- cc[i]*B[i]*(-1/(alpha)+(Z[i]^2)/(alpha))+
        (1-cc[i])*B[i]*(Z[i]*pdfZ[i])/(alpha*S_y[i])
      Aux4 <- matrix(c(Aux1,Aux2,Aux3), (p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)

    }

  }

  if(dist == "lognormal" && link=="cauchit"){
    mu <- (x %*% beta)
    theta<-(1/pi)*atan(w%*%eta)+0.5
    derv<-(1/pi)*(1/(1+(w%*%eta)^2))

    Z <- (log(y) - mu)/alpha
    pdfZ <- dnorm(Z)
    cdfZ <- pnorm(Z)
    S_y <- 1 - cdfZ

    MI<-matrix(0,(p+q+1),(p+q+1))
    for (i in 1:n){
      S_y[i] <- max(S_y[i], .Machine$double.eps)
      Aux1 <- (cc[i]*B[i]*(Z[i]/alpha)+
                 (1-cc[i])*B[i]*(pdfZ[i]/(S_y[i]*alpha)))%*%(x[i,])
      Aux2<- (B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])
      Aux3 <- cc[i]*B[i]*(-1/(alpha)+(Z[i]^2)/(alpha))+
        (1-cc[i])*B[i]*(Z[i]*pdfZ[i])/(alpha*S_y[i])
      Aux4 <- matrix(c(Aux1,Aux2,Aux3), (p+q+1),1)
      MI<- MI+Aux4%*%t(Aux4)

    }

  }

  if(dist == "loglogistic" && link=="probit"){
    lambda<- exp(x %*% beta)
    theta<-pnorm(w %*% eta)
    pdfN<-dnorm(w %*% eta)
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      y_alpha <- (y / lambda)^alpha  # Transformed survival time
      # Log-Logistic survival function
      S_y <- 1 / (1 + y_alpha)
      # Log-Logistic density function
      f_y <- (alpha / lambda) * (y / lambda)^(alpha - 1) / (1 + y_alpha)^2
      # Avoid division by zero issues
      S_y <- max(S_y, .Machine$double.eps)

      # Gradient components
      Aux1 <- (cc[i] * B[i] * (-alpha * (1 - y_alpha[i]) / (1 + y_alpha[i])) -
                 (1 - cc[i]) * B[i] * (alpha * y_alpha[i] / (1 + y_alpha[i]))) * x[i,]

      Aux2<- (B[i]*pdfN[i]/(theta[i])-(1-B[i])*pdfN[i]/(1-theta[i]))%*%(w[i,])

      Aux3 <- (cc[i] * B[i] * (1 / alpha + log(y[i]) - log(lambda[i]) -
                                 2* (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))) -
        ((1 - cc[i]) * B[i] * (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))

      # Combine into one vector
      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)

      # Update observed Fisher information matrix
      MI <- MI + Aux4 %*% t(Aux4)
    }
  }

  if(dist == "loglogistic" && link=="logit"){
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))

    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      y_alpha <- (y / lambda)^alpha  # Transformed survival time

      # Log-Logistic survival function
      S_y <- 1 / (1 + y_alpha)

      # Log-Logistic density function
      f_y <- (alpha / lambda) * (y / lambda)^(alpha - 1) / (1 + y_alpha)^2

      # Avoid division by zero issues
      S_y <- max(S_y, .Machine$double.eps)

      # Gradient components
      Aux1 <- (cc[i] * B[i] * (-alpha * (1 - y_alpha[i]) / (1 + y_alpha[i])) -
                 (1 - cc[i]) * B[i] * (alpha * y_alpha[i] / (1 + y_alpha[i]))) * x[i,]

      Aux2<- (B[i]*(1-theta[i])-(1-B[i])*theta[i])%*%(w[i,])

      Aux3 <- (cc[i] * B[i] * (1 / alpha + log(y[i]) - log(lambda[i]) -
                                 2* (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))) -
        ((1 - cc[i]) * B[i] * (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))

      # Combine into one vector
      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)

      # Update observed Fisher information matrix
      MI <- MI + Aux4 %*% t(Aux4)
    }
  }

  if(dist == "loglogistic" && link=="plogit"){
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))^tau
    derv<-tau*exp(-w %*% eta)/(1+exp(-w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      y_alpha <- (y / lambda[i])^alpha  # Transformed survival time

      # Log-Logistic survival function
      S_y <- 1 / (1 + y_alpha)

      # Log-Logistic density function
      f_y <- (alpha / lambda[i]) * (y / lambda[i])^(alpha - 1) / (1 + y_alpha)^2

      # Avoid division by zero issues
      S_y <- max(S_y, .Machine$double.eps)

      # Gradient components
      Aux1 <- (cc[i] * B[i] * (-alpha * (1 - y_alpha[i]) / (1 + y_alpha[i])) -
                 (1 - cc[i]) * B[i] * (alpha * y_alpha[i] / (1 + y_alpha[i]))) * x[i,]

      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])

      Aux3 <- (cc[i] * B[i] * (1 / alpha + log(y[i]) - log(lambda[i]) -
                                 2* (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))) -
        ((1 - cc[i]) * B[i] * (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))

      # Combine into one vector
      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)

      # Update observed Fisher information matrix
      MI <- MI + Aux4 %*% t(Aux4)
    }
  }

  if(dist == "loglogistic" && link=="rplogit"){
    lambda<- exp(x %*% beta)
    theta<-1-(1/(1+exp(w %*% eta)))^tau
    derv<-tau*exp(w %*% eta)/(1+exp(w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      y_alpha <- (y / lambda)^alpha  # Transformed survival time

      # Log-Logistic survival function
      S_y <- 1 / (1 + y_alpha)

      # Log-Logistic density function
      f_y <- (alpha / lambda) * (y / lambda)^(alpha - 1) / (1 + y_alpha)^2

      # Avoid division by zero issues
      S_y <- max(S_y, .Machine$double.eps)

      # Gradient components
      Aux1 <- (cc[i] * B[i] * (-alpha * (1 - y_alpha[i]) / (1 + y_alpha[i])) -
                 (1 - cc[i]) * B[i] * (alpha * y_alpha[i] / (1 + y_alpha[i]))) * x[i,]

      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])

      Aux3 <- (cc[i] * B[i] * (1 / alpha + log(y[i]) - log(lambda[i]) -
                                 2* (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))) -
        ((1 - cc[i]) * B[i] * (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))

      # Combine into one vector
      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)

      # Update observed Fisher information matrix
      MI <- MI + Aux4 %*% t(Aux4)
    }
  }

  if(dist == "loglogistic" && link=="cauchit"){
    lambda<- exp(x %*% beta)
    theta<-(1/pi)*atan(w%*%eta)+0.5
    derv<-(1/pi)*(1/(1+(w%*%eta)^2))
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      y_alpha <- (y / lambda)^alpha  # Transformed survival time

      # Log-Logistic survival function
      S_y <- 1 / (1 + y_alpha)

      # Log-Logistic density function
      f_y <- (alpha / lambda) * (y / lambda)^(alpha - 1) / (1 + y_alpha)^2

      # Avoid division by zero issues
      S_y <- max(S_y, .Machine$double.eps)

      # Gradient components
      Aux1 <- (cc[i] * B[i] * (-alpha * (1 - y_alpha[i]) / (1 + y_alpha[i])) -
                 (1 - cc[i]) * B[i] * (alpha * y_alpha[i] / (1 + y_alpha[i]))) * x[i,]

      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])

      Aux3 <- (cc[i] * B[i] * (1 / alpha + log(y[i]) - log(lambda[i]) -
                                 2* (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))) -
        ((1 - cc[i]) * B[i] * (y_alpha[i] * log(y[i] / lambda[i])) / (1 + y_alpha[i]))

      # Combine into one vector
      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)

      # Update observed Fisher information matrix
      MI <- MI + Aux4 %*% t(Aux4)
    }
  }

  if(dist == "invgauss" && link == "logit"){
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))

    if(length(which(B == 0)) > 0) B[which(B == 0)] <- .Machine$double.xmin

    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      ##--- (1) base inverse‐Gaussian score pieces
      dlogf_alpha <- 1/(2*alpha)-
        (y[i] - lambda[i])^2 / (2 * lambda[i]^2 * y[i])

      dlogf_beta  <- alpha*((y[i] - lambda[i])/(lambda[i]^2))

      ##--- (2) compute the z’s
      z1 <- sqrt(alpha/y[i])*(y[i]/lambda[i] - 1)
      z2 <- sqrt(alpha/y[i])*(y[i]/lambda[i] + 1)

      SS <- actuar::pinvgauss(q = y[i], mean = lambda[i],
                              shape = alpha,
                              lower.tail = FALSE)
      SS <- max(SS, .Machine$double.xmin)

      z1pdf <- max(dnorm(z1), .Machine$double.xmin)
      z2pdf <- max(dnorm(-z2), .Machine$double.xmin)
      z1cdf <- max(pnorm(z1), .Machine$double.xmin)
      z2cdf <- max(pnorm(-z2), .Machine$double.xmin)


      dS_alpha <- -z1pdf * (y[i]/lambda[i] - 1)/(2*sqrt(alpha*y[i])) -
        exp(2*alpha/lambda[i])*(2/lambda[i])*z2cdf +
        exp(2*alpha/lambda[i])*z2pdf*(y[i]/lambda[i] + 1)/(2*sqrt(alpha*y[i]))

      dlogS_alpha <- dS_alpha/SS


      dS_lambda <- z1pdf*(sqrt(alpha*y[i]))/(lambda[i]^2) +
        exp(2*alpha/lambda[i])*2*alpha/(lambda[i]^2)*z2cdf -
        exp(2*alpha/lambda[i]) * z2pdf *(sqrt(alpha * y[i]))/(lambda[i]^2)

      dS_beta <- dS_lambda * lambda[i]
      dlogS_beta <- dS_beta/SS


      ##--- (8) mixture‐model contributions
      Aux1 <- (cc[i]*B[i]*dlogf_beta+
                 (1 - cc[i])*B[i]*dlogS_beta)*x[i,]

      Aux2<- (B[i]*(1-theta[i])-(1-B[i])*theta[i])%*%(w[i,])

      Aux3 <- cc[i] * B[i] * dlogf_alpha +
        (1 - cc[i]) * B[i] * dlogS_alpha

      ##--- (9) update Fisher information
      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)
      if(sum(is.na(Aux4)) > 0) Aux4 <- matrix(0, (p + q + 1), 1)
      MI   <- MI + Aux4 %*% t(Aux4)
    }
  }


  if(dist == "invgauss" && link == "probit"){
    lambda<- exp(x %*% beta)
    theta<-pnorm(w %*% eta)
    pdfN<-dnorm(w %*% eta)
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n) {
      ##--- (1) base inverse‐Gaussian score pieces
      dlogf_alpha <- 1/(2*alpha) -
        (y[i] - lambda[i])^2 / (2 * lambda[i]^2 * y[i])
      dlogf_beta  <- (y[i] - lambda[i]) / (lambda[i]^2) * alpha

      ##--- (2) compute the z’s
      z1 <- sqrt(alpha/y[i])*(y[i]/lambda[i] - 1)
      z2 <- sqrt(alpha/y[i])*(y[i]/lambda[i] + 1)

      SS <- actuar::pinvgauss(q = y[i], mean = lambda[i], shape = alpha, lower.tail = FALSE)
      if(SS == 0) SS <- .Machine$double.xmin

      z1pdf <- dnorm(z1)
      z2pdf <- dnorm(-z2)
      z1cdf <- pnorm(z1)
      z2cdf <- pnorm(-z2)

      dS_alpha <- -z1pdf * (y[i]/lambda[i] - 1)/(2*sqrt(alpha*y[i])) -
        exp(2*alpha/lambda[i])*(2/lambda[i])*z2cdf +
        exp(2*alpha/lambda[i])*z2pdf*(y[i]/lambda[i] + 1)/(2*sqrt(alpha * y[i]))

      dlogS_alpha <- dS_alpha/SS


      dS_lambda <- z1pdf*(sqrt(alpha*y[i]))/(lambda[i]^2) +
        exp(2*alpha/lambda[i])*2*alpha/(lambda[i]^2)*z2cdf -
        exp(2*alpha/lambda[i]) * z2pdf *(sqrt(alpha * y[i]))/(lambda[i]^2)

      dS_beta <- dS_lambda * lambda[i]
      dlogS_beta <- dS_beta/SS

      ##--- (8) mixture‐model contributions
      Aux1 <- ( cc[i] * B[i] * dlogf_beta
                + (1 - cc[i]) * B[i] * dlogS_beta
      ) * x[i,]

      Aux2 <- ( B[i]*pdfN[i]/theta[i]
                - (1-B[i])*pdfN[i]/(1-theta[i])
      ) %*% w[i,]

      Aux3 <- cc[i] * B[i] * dlogf_alpha +
        (1 - cc[i]) * B[i] * dlogS_alpha

      ##--- (9) update Fisher information
      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)
      if(sum(is.na(Aux4)) > 0) Aux4 <- matrix(0, (p + q + 1), 1)
      MI   <- MI + Aux4 %*% t(Aux4)
    }
  }

  if(dist == "invgauss" && link == "plogit"){
    lambda<- exp(x %*% beta)
    theta<-1/(1+exp(-w %*% eta))^tau
    derv<-tau*exp(-w %*% eta)/(1+exp(-w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      ##--- (1) base inverse‐Gaussian score pieces
      dlogf_alpha <- 1/(2*alpha) -
        (y[i] - lambda[i])^2 / (2 * lambda[i]^2 * y[i])
      dlogf_beta  <- (y[i] - lambda[i]) / (lambda[i]^2) * alpha

      ##--- (2) compute the z’s
      z1 <- sqrt(alpha/y[i])*(y[i]/lambda[i] - 1)
      z2 <- sqrt(alpha/y[i])*(y[i]/lambda[i] + 1)

      SS <- actuar::pinvgauss(q = y[i], mean = lambda[i], shape = alpha, lower.tail = FALSE)
      if(SS == 0) SS <- .Machine$double.xmin

      z1pdf <- dnorm(z1)
      z2pdf <- dnorm(-z2)
      z1cdf <- pnorm(z1)
      z2cdf <- pnorm(-z2)

      dS_alpha <- -z1pdf * (y[i]/lambda[i] - 1)/(2*sqrt(alpha*y[i])) -
        exp(2*alpha/lambda[i])*(2/lambda[i])*z2cdf +
        exp(2*alpha/lambda[i])*z2pdf*(y[i]/lambda[i] + 1)/(2*sqrt(alpha * y[i]))

      dlogS_alpha <- dS_alpha/SS


      dS_lambda <- z1pdf*(sqrt(alpha*y[i]))/(lambda[i]^2) +
        exp(2*alpha/lambda[i])*2*alpha/(lambda[i]^2)*z2cdf -
        exp(2*alpha/lambda[i]) * z2pdf *(sqrt(alpha * y[i]))/(lambda[i]^2)

      dS_beta <- dS_lambda * lambda[i]
      dlogS_beta <- dS_beta/SS

      ##--- (8) mixture‐model contributions
      Aux1 <- ( cc[i] * B[i] * dlogf_beta
                + (1 - cc[i]) * B[i] * dlogS_beta
      ) * x[i,]

      Aux2 <- (B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])

      Aux3 <- cc[i] * B[i] * dlogf_alpha +
        (1 - cc[i]) * B[i] * dlogS_alpha

      ##--- (9) update Fisher information
      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)
      if(sum(is.na(Aux4)) > 0) Aux4 <- matrix(0, (p + q + 1), 1)
      MI   <- MI + Aux4 %*% t(Aux4)
    }
  }

  if(dist == "invgauss" && link == "rplogit"){
    lambda<- exp(x %*% beta)
    theta<-1-(1/(1+exp(w %*% eta)))^tau
    derv<-tau*exp(w %*% eta)/(1+exp(w %*% eta))^(tau+1)
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      ##--- (1) base inverse‐Gaussian score pieces
      dlogf_alpha <- 1/(2*alpha) -
        (y[i] - lambda[i])^2 / (2 * lambda[i]^2 * y[i])
      dlogf_beta  <- (y[i] - lambda[i]) / (lambda[i]^2) * alpha

      ##--- (2) compute the z’s
      z1 <- sqrt(alpha/y[i])*(y[i]/lambda[i] - 1)
      z2 <- sqrt(alpha/y[i])*(y[i]/lambda[i] + 1)

      SS <- actuar::pinvgauss(q = y[i], mean = lambda[i], shape = alpha, lower.tail = FALSE)
      if(SS == 0) SS <- .Machine$double.xmin

      z1pdf <- dnorm(z1)
      z2pdf <- dnorm(-z2)
      z1cdf <- pnorm(z1)
      z2cdf <- pnorm(-z2)

      dS_alpha <- -z1pdf * (y[i]/lambda[i] - 1)/(2*sqrt(alpha*y[i])) -
        exp(2*alpha/lambda[i])*(2/lambda[i])*z2cdf +
        exp(2*alpha/lambda[i])*z2pdf*(y[i]/lambda[i] + 1)/(2*sqrt(alpha * y[i]))

      dlogS_alpha <- dS_alpha/SS


      dS_lambda <- z1pdf*(sqrt(alpha*y[i]))/(lambda[i]^2) +
        exp(2*alpha/lambda[i])*2*alpha/(lambda[i]^2)*z2cdf -
        exp(2*alpha/lambda[i]) * z2pdf *(sqrt(alpha * y[i]))/(lambda[i]^2)

      dS_beta <- dS_lambda * lambda[i]
      dlogS_beta <- dS_beta/SS

      ##--- (8) mixture‐model contributions
      Aux1 <- ( cc[i] * B[i] * dlogf_beta
                + (1 - cc[i]) * B[i] * dlogS_beta
      ) * x[i,]

      Aux2 <- (B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])

      Aux3 <- cc[i] * B[i] * dlogf_alpha +
        (1 - cc[i]) * B[i] * dlogS_alpha

      ##--- (9) update Fisher information
      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)
      if(sum(is.na(Aux4)) > 0) Aux4 <- matrix(0, (p + q + 1), 1)
      MI   <- MI + Aux4 %*% t(Aux4)
    }
  }

  if(dist == "invgauss" && link == "cauchit"){
    lambda<- exp(x %*% beta)
    theta<-(1/pi)*atan(w%*%eta)+0.5
    derv<-(1/pi)*(1/(1+(w%*%eta)^2))
    MI<-matrix(0,(p+q+1),(p+q+1))

    for (i in 1:n){
      ##--- (1) base inverse‐Gaussian score pieces
      dlogf_alpha <- 1/(2*alpha) -
        (y[i] - lambda[i])^2 / (2 * lambda[i]^2 * y[i])
      dlogf_beta  <- (y[i] - lambda[i]) / (lambda[i]^2) * alpha

      ##--- (2) compute the z’s
      z1 <- sqrt(alpha/y[i])*(y[i]/lambda[i] - 1)
      z2 <- sqrt(alpha/y[i])*(y[i]/lambda[i] + 1)

      SS <- actuar::pinvgauss(q = y[i], mean = lambda[i], shape = alpha, lower.tail = FALSE)
      if(SS == 0) SS <- .Machine$double.xmin

      z1pdf <- dnorm(z1)
      z2pdf <- dnorm(-z2)
      z1cdf <- pnorm(z1)
      z2cdf <- pnorm(-z2)

      dS_alpha <- -z1pdf * (y[i]/lambda[i] - 1)/(2*sqrt(alpha*y[i])) -
        exp(2*alpha/lambda[i])*(2/lambda[i])*z2cdf +
        exp(2*alpha/lambda[i])*z2pdf*(y[i]/lambda[i] + 1)/(2*sqrt(alpha * y[i]))

      dlogS_alpha <- dS_alpha/SS


      dS_lambda <- z1pdf*(sqrt(alpha*y[i]))/(lambda[i]^2) +
        exp(2*alpha/lambda[i])*2*alpha/(lambda[i]^2)*z2cdf -
        exp(2*alpha/lambda[i]) * z2pdf *(sqrt(alpha * y[i]))/(lambda[i]^2)

      dS_beta <- dS_lambda * lambda[i]
      dlogS_beta <- dS_beta/SS

      Aux1 <- ( cc[i] * B[i] * dlogf_beta
                + (1 - cc[i]) * B[i] * dlogS_beta
      ) * x[i,]

      Aux2<-(B[i]*derv[i]/(theta[i])-(1-B[i])*derv[i]/(1-theta[i]))%*%(w[i,])

      Aux3 <- cc[i] * B[i] * dlogf_alpha +
        (1 - cc[i]) * B[i] * dlogS_alpha

      Aux4 <- matrix(c(Aux1, Aux2, Aux3), (p + q + 1), 1)
      if(sum(is.na(Aux4)) > 0) Aux4 <- matrix(0, (p + q + 1), 1)
      MI   <- MI + Aux4 %*% t(Aux4)
    }
  }

  ep.aux <- sqrt(diag(solve(MI)))
  ep.beta <- ep.aux[1:p]
  ep.eta <- ep.aux[(p+1):(p+q)]

  if(dist %in% c("exponential", "rayleigh")){
    ep.final <- c(ep.beta, ep.eta)
  }else{
    ep.alpha <- ep.aux[(p+q+1)]
    ep.final <- c(ep.alpha, ep.beta, ep.eta)
  }

  return(ep.final)
}
#
format_with_dash <- function(x){
  if (is.character(x)){
    x_numeric <- as.numeric(ifelse(x == "-", NA, x))
    formatted <- ifelse(is.na(x_numeric), "-", format(round(x_numeric, 3),
                                                      nsmall = 3))
    return(formatted)
  } else {
    return(round(x,3))
  }
}
#' Fit a Mixture Cure Rate (MCR) Survival Model
#'
#' Fits a cure rate model using a flexible link function and a variety of survival distributions. The model accounts for a cured fraction through a logistic-type link and estimates the model via an EM-like algorithm.
#'
#' @import Formula
#' @importFrom survival Surv survreg survfit
#' @import knitr
#' @import flexsurv
#' @import tibble
#' @import stats
#' @importFrom actuar dinvgauss pinvgauss
#'
#' @param formula A two-part formula of the form \code{Surv(time, status) ~ x | w}, where \code{x} are covariates for the survival part, and \code{w} are covariates for the cure fraction.
#' @param data A data frame containing the variables in the model.
#' @param dist A character string indicating the baseline distribution. Supported values are \code{"weibull"}, \code{"exponential"}, \code{"rayleigh"}, \code{"lognormal"}, \code{"loglogistic"}, and \code{"invgauss"}.
#' @param link A character string specifying the link function for the cure fraction. Options are \code{"logit"}, \code{"probit"}, \code{"plogit"}, \code{"rplogit"}, and \code{"cauchit"}.
#' @param tau A numeric value used when \code{link = "plogit"} or \code{"rplogit"}. Defaults to 1.
#' @param maxit Maximum number of iterations for the EM-like algorithm. Defaults to 1000.
#' @param tol Convergence tolerance. Defaults to 1e-5.
#'
#' @return An object of class \code{"MCR"}, which is a list containing:
#' \item{coefficients}{Estimated regression coefficients for the survival part.}
#' \item{coefficients_cure}{Estimated coefficients for the cure part.}
#' \item{scale}{Estimated scale parameter of the baseline distribution.}
#' \item{loglik}{Final log-likelihood value.}
#' \item{n}{Number of observations used in the model.}
#' \item{deleted}{Number of incomplete cases removed before fitting.}
#' \item{ep}{Estimated standard errors.}
#' \item{iter}{Number of iterations used for convergence.}
#' \item{dist}{Distribution used.}
#' \item{link}{Link function used.}
#' \item{tau}{Tau parameter used (if applicable).}
#'
#' @examples
#' require(EMGCR)
#'
#' data(liver2)
#' names(liver2)
#' liver2$sex <- factor(liver2$sex)
#' liver2$grade <- factor(liver2$grade)
#' liver2$radio <- factor(liver2$radio)
#' liver2$chemo <- factor(liver2$chemo)
#' str(liver2)
#' model <- MCRfit(
#'   survival::Surv(time, status) ~ age + sex + grade + radio + chemo |
#'     age + medh + grade + radio + chemo,
#'   dist = "loglogistic",
#'   link = "plogit",
#'   tau = 0.15,
#'   data = liver2
#' )
#' summary(model)
#'
#' @export
#'
MCRfit<-function(formula,data,dist="weibull",
                link="logit",tau=1,
                maxit = 1E3, tol = 1E-5){


  #

  mf1 <- model.frame(Formula(formula), data = data)
  x <- model.matrix(Formula(formula), data = mf1, rhs = 1)
  w <- model.matrix(Formula(formula), data = mf1, rhs = 2)

  #

  mf <- model.frame(Formula(formula), data = data)
  model.aux <- model.response(mf)
  cc <- model.aux[, "status"]
  y <- model.aux[, "time"]


  if(dist != "invgauss"){
  fit0 <- survreg(survival::Surv(y, cc) ~ x - 1, dist = dist)
  beta0 <- as.vector(fit0$coeff)
  eta0 <- as.vector(-glm(cc ~ w - 1, family = binomial)$coeff)

  if (dist %in% c("exponential", "rayleigh")) {
    alpha0 <- 1/fit0$scale
    beta0 <- beta0*alpha0
    para1 <- c(beta0, eta0)
  } else if (dist == "weibull") {
    alpha0 <- 1/fit0$scale
    beta0 <- beta0*alpha0
    para1 <- c(alpha0, beta0, eta0)
  } else if (dist == "lognormal") {
    alpha0 <- fit0$scale
    para1 <- c(alpha0, beta0, eta0)
  } else if (dist == "loglogistic") {
    alpha0 <- 1/fit0$scale
    para1 <- c(alpha0, beta0, eta0)
  } else {
    stop("Unsupported distribution function.")
  }

  #
  lambda <- exp(x%*%beta0)
  if(link=="logit"){
    theta<-1/(1+exp(-w%*%eta0))
  }else if(link=="probit"){
      theta<-pnorm(w%*%eta0)
  } else  if(link=="plogit"){
        theta<-(1/(1+exp(-w%*%eta0)))^tau
  } else if(link=="rplogit"){
    theta<-1-(1/(1+exp(w %*% eta0)))^tau
  } else if(link=="cauchit"){
    theta<-(1/pi)*atan(w%*%eta0)+0.5
  } else{
    stop("Unsupported link function.")
        }

  iter <- 0
  criteria <- 1

  while((criteria > tol) && (iter <= maxit)){

    if(dist=="exponential"){
      aux<- pexp(q=y,rate=lambda,lower.tail = FALSE) #para1[1]
    } else if(dist=="rayleigh"){
      aux<- pweibull(q=y,shape=2,scale=lambda**(-1/2),lower.tail = FALSE) #para1[1]
    } else if(dist=="weibull"){
      aux<- pweibull(q=y,shape=para1[1],scale=(lambda)**(-1/para1[1]),
                     lower.tail = FALSE) #para1[1]
    } else if(dist=="lognormal"){
      aux<- plnorm(q=y,meanlog = -log(lambda),
                   sdlog = para1[1],lower.tail = FALSE) #para1[1]
    } else if(dist == "loglogistic"){
      aux <- flexsurv::pllogis(q=y, shape=para1[1],scale=lambda,
                               lower.tail = FALSE) #para1[1]
    }


    B <- cc+(1-cc)*theta*aux/(1-theta+theta*aux)


    if(link=="logit"){eta<- suppressWarnings(as.vector(glm(c(B)~w-1,
                        family=binomial(link=logit))$coeff))
    } else if(link=="probit") {
      eta<- suppressWarnings(as.vector(glm(c(B)~w-1,
                        family=binomial(link=probit))$coeff))
      } else if(link=="plogit") {
        eta<- suppressWarnings(as.vector(glm(c(B)~w-1,
                        family=binomial(link =plogit(tau = tau)),
                        epsilon = 1E-8, maxit = 1E3)$coeff))
      } else if(link=="cauchit"){
        eta<- suppressWarnings(as.vector(glm(c(B)~w-1,
                        family=binomial(link=cauchit))$coeff))
      } else if(link=="rplogit"){
        eta<- suppressWarnings(as.vector(glm(c(B)~w-1,
                        family=binomial(link =rplogit(tau = tau)),
                        epsilon = 1E-8, maxit = 1E3)$coeff))
        }

    if(length(which(B == 0)) > 0) B[which(B == 0)] <- .Machine$double.xmin

    fit <- survreg(survival::Surv(y, cc) ~ x - 1,
                             weights = c(B), dist = dist)

    if (dist %in% c("exponential", "rayleigh")) {
      alpha <- 1/fit$scale
      beta <- -fit$coeff*alpha
    } else if (dist == "weibull") {
      alpha <- 1/fit$scale
      beta <- -fit$coeff*alpha
    } else if (dist == "lognormal") {
      alpha <- fit$scale
      beta <- -fit$coeff
    } else if (dist == "loglogistic") {
      alpha <- 1/fit$scale
      beta <- fit$coeff
    } else {
      stop("Unsupported distribution.")
    }

    lambda <- exp(x %*% beta)
    if(link=="logit"){
      theta<-1/(1+exp(-w %*% eta))
    }else if(link=="probit") {
        theta<-pnorm(w%*%eta)
    }else if(link=="plogit"){
          theta<-(1/(1+exp(-w %*% eta)))^tau
    } else if(link =="rplogit"){
      theta<-1-(1/(1+exp(w %*% eta)))^tau
    } else if(link =="cauchit"){
      theta<-(1/pi)*atan(w%*%eta)+0.5
          }

    if (dist %in% c("exponential", "rayleigh")){
      para2 <- c(beta, eta)
    }else{
      para2 <- c(alpha, beta, eta)
    }

    criteria <- (para2-para1)%*%(para2-para1)

    iter <- iter+1
    para1 <- para2
  }
  }else{
    pp <- ncol(x)
    fit.ini <- survreg(survival::Surv(y,cc) ~ x - 1, dist = "weibull")
    alpha <- 1/fit.ini$scale
    beta <- -fit.ini$coeff*alpha
    names(beta) <- NULL
    eta <- as.vector(glm(cc~w-1,family=binomial)$coeff)

    para1 <- c(alpha, beta, eta)

    #
    lambda<- exp(x%*%beta)
    if(link=="logit"){
      theta<-1/(1+exp(-w%*%eta))
    }else if(link=="probit"){
        theta<-pnorm(w%*%eta)
    }else if(link=="plogit"){
          theta<-(1/(1+exp(-w%*%eta)))^tau
    }else if(link=="rplogit"){
            theta<-1-(1/(1+exp(w %*% eta)))^tau
    }else if(link=="cauchit"){
              theta<-(1/pi)*atan(w%*%eta)+0.5
    }

    iter <- 0
    criteria <- 1

    while((criteria > tol) && (iter <= maxit)){

      aux<- actuar::pinvgauss(q = y, mean=lambda,shape = para1[1],
                              lower.tail = FALSE,log.p = FALSE)

      B <- cc+(1-cc)*theta*aux/(1-theta+theta*aux)
      if(link=="logit"){
        eta<- suppressWarnings(as.vector(glm(c(B)~ w-1,
                                             family=binomial)$coeff))
        }else if(link=="probit"){
          eta<- suppressWarnings(as.vector(glm(c(B)~ w-1,
                                               family=binomial(link=probit))$coeff))
          }else if(link=="plogit"){
            eta<- suppressWarnings(as.vector(glm(c(B)~ w-1,
                                                 family=binomial(link =plogit(tau = tau)),epsilon = 1e-08, maxit = 1000)$coeff))
            }else if(link=="cauchit"){
              eta<- suppressWarnings(as.vector(glm(c(B)~ w-1,family=binomial(link=cauchit))$coeff))
              }else if(link=="rplogit"){
                eta<- suppressWarnings(as.vector(glm(c(B)~ w-1,family=binomial(link =rplogit(tau = tau)), epsilon = 1e-08, maxit = 1000)$coeff))
                }else if(length(which(B == 0)) > 0){
                  B[which(B == 0)] <- .Machine$double.xmin
                  }


      #grad <- gradInvGauss(params=c(beta,alpha), time=y, status=cc,
      #                  x=x,B=B)
      #print(grad)

      fit<- optim(par=c(beta,alpha),
                      method = "L-BFGS-B", fn=likeInvGauss,
                      lower = c(rep(-Inf, pp),1e-06),
                      upper = c(rep(Inf, pp),Inf),
                      time=y, status=cc, x=x, B=B,
                      hessian=TRUE,
                      control = list(maxit = 30000,
                                     temp = 2000, trace = FALSE,
                                     REPORT = 500))

        #print(fit$par)

        beta<- fit$par[1:pp]
        alpha<- fit$par[pp+1]

        lambda <- exp(x %*% beta)
        if(link=="logit"){
          theta<-1/(1+exp(-w %*% eta))
        }else if(link=="probit"){
          theta<-pnorm(w%*%eta)
        }else if(link=="plogit"){
          theta<-(1/(1+exp(-w %*% eta)))^tau
        }else if(link=="rplogit"){
          theta<-1-(1/(1+exp(w %*% eta)))^tau
        }else if(link=="cauchit"){
          theta<-(1/pi)*atan(w%*%eta)+0.5
        }

      para2 <- c(alpha,beta,eta)

      iter <- iter+1

      #print(paste("Iteration :", iter,sep=""))

      criteria <- (para2-para1)%*%(para2-para1)

      para1 <- para2
    }
  }


  loglink <- likeMRC(y=y,cc=cc,x=x,w=w,alpha=alpha,beta=beta,
                      eta=eta,dist=dist,link=link,tau=tau)

  ep <- epMCR(y=y,cc=cc,x=x,w=w,B=B,alpha=alpha,beta=beta,
              eta=eta,tau=tau,dist=dist,link=link)

  names(beta) <- colnames(x)
  names(eta) <- colnames(w)

  fit.MCR <- list(
    call = match.call(),
    formula = formula,
    coefficients = beta,
    coefficients_cure = eta,
    scale = alpha,
    loglik = loglink,
    n = nrow(data),
    deleted = sum(!complete.cases(model.frame(Formula(formula), data = data))),
    ep = ep,
    iter = iter,
    dist = dist,
    link = link,
    tau = tau,
    data = data
  )
  fit.MCR$data <- data
  class(fit.MCR) <- "MCR"
  return(fit.MCR)

}
#' @export
print.MCR <- function(x, ...) {
  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients (survival part):\n")
  print(round(x$coefficients, 4))

  cat("\nCoefficients (cure part):\n")
  print(round(x$coefficients_cure, 4))

  cat("\nScale:\n")
  print(round(x$scale, 4))

  cat("\nLog-likelihood:", round(x$loglik, 4), "\n")
}

#' @export
summary.MCR <- function(object, ...) {
  y <- model.response(model.frame(Formula(object$formula), data = eval(object$call$data)))
  cc <- y[, "status"]
  y <- y[, "time"]

  x <- model.matrix(Formula(object$formula), data = eval(object$call$data), rhs = 1)
  w <- model.matrix(Formula(object$formula), data = eval(object$call$data), rhs = 2)

  coef_s <- object$coefficients
  coef_cure <- object$coefficients_cure
  scale <- object$scale
  ep <- object$ep
  iter <- object$iter
  value <- object$loglik
  dist <- object$dist
  tau <- object$tau
  n <- object$n

  n_alpha <- if (dist %in% c("weibull", "lognormal", "loglogistic","invgauss")) 1 else 0
  std_alpha <- if (n_alpha == 1) ep[1] else NULL
  std_beta <- ep[(1 + n_alpha):(n_alpha + length(coef_s))]
  std_eta <- ep[(n_alpha + length(coef_s) + 1):length(ep)]


  z_beta <- coef_s / std_beta
  p_beta <- 2 * (1 - pnorm(abs(z_beta)))

  z_eta <- coef_cure / std_eta
  p_eta <- 2 * (1 - pnorm(abs(z_eta)))

  coef_surv <- data.frame(
    Value = coef_s,
    `Std. Error` = std_beta,
    z = z_beta,
    p = format.pval(p_beta, digits = 2, eps = 2e-16),
    row.names = colnames(x),
    check.names = FALSE
  )

  coef_cure <- data.frame(
    Value = coef_cure,
    `Std. Error` = std_eta,
    z = z_eta,
    p = format.pval(p_eta, digits = 2 , eps = 2e-16),
    row.names = colnames(w),
    check.names = FALSE
  )

  if (!is.null(object$scale) && length(object$scale) == 1 && object$scale != 1 && object$scale != 2) {
    std_alpha <- object$ep[1]
    coef_scale <- data.frame(
      Estimate = object$scale,
      Std.Error = std_alpha,
      row.names = "alpha"
    )
  } else {
    coef_scale <- NULL
  }

  out <- list(
    call = object$call,
    dist = dist,
    loglik = value,
    AIC = -2 * value + 2 * length(ep),
    BIC = -2 * value + log(n) * length(ep),
    tau = tau,
    iter = iter,
    coef_surv = coef_surv,
    coef_cure = coef_cure,
    scale = coef_scale
  )

  class(out) <- "summary.MCR"
  return(out)
}
#' @export
print.summary.MCR <- function(x, digits = 5, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nDistribution:", x$dist, "\n")
  cat("Log-Likelihood:", formatC(x$loglik, digits = digits, format = "f"), "\n")
  cat("AIC:", formatC(x$AIC, digits = digits, format = "f"), "\n")
  cat("BIC:", formatC(x$BIC, digits = digits, format = "f"), "\n")
  cat("tau:", x$tau, "\n")
  cat("Number of Iterations:", x$iter, "\n\n")

  cat("Coefficients (Survival part):\n")
  print(format(x$coef_surv, digits = digits, nsmall = digits), quote = FALSE)

  cat("\nCoefficients (Uncured part):\n")
  print(format(x$coef_cure, digits = digits, nsmall = digits), quote = FALSE)

  if (!is.null(x$scale)) {
    cat("\nScale:\n")
    print(format(x$scale, digits = digits, nsmall = digits), quote = FALSE)
  }
}
