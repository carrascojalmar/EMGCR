# Exponential distribution
exponential_pdf <- function(x, alpha=1, lambda) {
  pdf <- alpha*lambda*x^(alpha-1)*exp(-lambda*x^alpha)
  return(pdf)
}

#
exponential_sf <- function(x, alpha=1, lambda) {
  sf <- exp(-lambda*x^alpha)
  return(sf)
}

# Rayleigh distribution
rayleigh_pdf <- function(x, alpha=2, lambda) {
  pdf <- alpha*lambda*x^(alpha-1)*exp(-lambda*x^alpha)
  return(pdf)
}

#
rayleigh_sf <- function(x, alpha=2, lambda) {
  sf <- exp(-lambda*x^alpha)
  return(sf)
}

# Weibull distribution
weibull_pdf <- function(x, alpha, lambda) {
  pdf <- alpha*lambda*x^(alpha-1)*exp(-lambda*x^alpha)
  return(pdf)
}

#
weibull_sf <- function(x, alpha, lambda) {
  sf <- exp(-lambda*x^alpha)
  return(sf)
}

# Log-normal distribution

lognormal_sf <- function(x,  alpha, mu) {
  sf <- 1-pnorm((log(x)-mu)/alpha)
  return(sf)
}

lognormal_pdf <- function(x, alpha, mu) {
  pdf <- dnorm(log(x), mu, alpha)/x
  return(pdf)
}

# Log-logistic distribution
loglogistic_sf <- function(x, alpha, lambda){
  sf <- 1/(1+(x/lambda)^alpha)
  return(sf)
}

# Log-logistic distribution
loglogistic_pdf <- function(x, alpha, lambda){
  pdf <- (alpha/lambda)*(x/lambda)^(alpha - 1)/(1+(x/lambda)^alpha)^2
  return(pdf)
}

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
  link <- list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, valideta = valideta, name = paste("plogit(tau=", tau, ")", sep=""))
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

likeMRC <- function(y, cc, x, w, alpha, beta, eta,
                    link=link, dist=dist, tau=tau){
  time <- y
  status <- cc

  if(dist=="exponential"){
    lambda <- exp(x%*%beta)
    if(link=="logit"){
      theta <- 1/(1+exp(-w%*%eta))
    }else if(link=="probit"){
      theta<-pnorm(w%*%eta)
    }else if(link=="plogit"){
      theta<-(1/(1+exp(-w%*%eta)))^tau
    }else if(link=="rplogit"){
      theta<-1-(1/(1+exp(w %*% eta)))^tau
    }else if(link=="cauchit"){
      theta<-(1/pi)*atan(w%*%eta)+0.5
    }

    ll <- numeric()
    for(i in 1:length(time)){
      aux1 <- (theta[i]*exponential_pdf(x=time[i],alpha=1,
                                    lambda=lambda[i]))^(status[i])
      aux2 <- (1-theta[i]+theta[i]*exponential_sf(x=time[i],alpha=1,
                                              lambda=lambda[i]))^(1-status[i])
      ll[i] <- aux1*aux2
    }
  }

  if(dist=="rayleigh"){
    lambda <- exp(x%*%beta)
    if(link=="logit"){
      theta <- 1/(1+exp(-w%*%eta))
    }else if(link=="probit"){
      theta<-pnorm(w%*%eta)
    }else if(link=="plogit"){
      theta<-(1/(1+exp(-w%*%eta)))^tau
    }else if(link=="rplogit"){
      theta<-1-(1/(1+exp(w %*% eta)))^tau
    }else if(link=="cauchit"){
      theta<-(1/pi)*atan(w%*%eta)+0.5
    }

    ll <- numeric()
    for(i in 1:length(time)){
      aux1 <- (theta[i]*rayleigh_pdf(x=time[i],alpha=2,
                                        lambda=lambda[i]))^(status[i])
      aux2 <- (1-theta[i]+theta[i]*rayleigh_sf(x=time[i],alpha=2,
                                                  lambda=lambda[i]))^(1-status[i])
      ll[i] <- aux1*aux2
    }
  }

  if(dist=="weibull"){
    lambda <- exp(x%*%beta)
    if(link=="logit"){
      theta <- 1/(1+exp(-w%*%eta))
    }else if(link=="probit"){
      theta<-pnorm(w%*%eta)
    }else if(link=="plogit"){
      theta<-(1/(1+exp(-w%*%eta)))^tau
    }else if(link=="rplogit"){
        theta<-1-(1/(1+exp(w %*% eta)))^tau
    }else if(link=="cauchit"){
          theta<-(1/pi)*atan(w%*%eta)+0.5
          }

    ll <- numeric()
    for(i in 1:length(time)){
      aux1 <- (theta[i]*weibull_pdf(x=time[i],alpha=alpha,
                                    lambda=lambda[i]))^(status[i])
      aux2 <- (1-theta[i]+theta[i]*weibull_sf(x=time[i],alpha=alpha,
                                              lambda=lambda[i]))^(1-status[i])
      ll[i] <- aux1*aux2
    }
  }

  if(dist=="lognormal"){
    lambda <- exp(x%*%beta)
    if(link=="logit"){
      theta <- 1/(1+exp(-w%*%eta))
    }else if(link=="probit"){
        theta<-pnorm(w%*%eta)
    }else if(link=="plogit"){
          theta<-(1/(1+exp(-w%*%eta)))^tau
    }else if(link=="rplogit"){
            theta<-1-(1/(1+exp(w %*% eta)))^tau
    } else if(link=="cauchit"){
              theta<-(1/pi)*atan(w%*%eta)+0.5
              }
    ll <- numeric()
    for(i in 1:length(time)){
      aux1 <- (theta[i]*lognormal_pdf(x=time[i],alpha=alpha,
              mu=-log(lambda[i])))^(status[i])
      aux2 <- (1-theta[i]+theta[i]*lognormal_sf(x=time[i],
              alpha=alpha,mu=-log(lambda[i])))^(1-status[i])
      ll[i] <- aux1*aux2
    }
  }

  if(dist=="loglogistic"){
    lambda <- exp(x%*%beta)
    if(link=="logit"){
      theta <- 1/(1+exp(-w%*%eta))
    }else if(link=="probit"){
      theta<-pnorm(w%*%eta)
    }else if(link=="plogit"){
      theta<-(1/(1+exp(-w%*%eta)))^tau
    }else if(link=="rplogit"){
      theta<-1-(1/(1+exp(w %*% eta)))^tau
    }else if(link=="cauchit"){
      theta<-(1/pi)*atan(w%*%eta)+0.5
    }

    ll <- numeric()
    for(i in 1:length(time)){
      aux1 <- (theta[i]*loglogistic_pdf(x=time[i],alpha=alpha,
                                        lambda=lambda[i]))^(status[i])
      aux2 <- (1-theta[i]+theta[i]*loglogistic_sf(x=time[i],alpha=alpha,
                                                  lambda=lambda[i]))^(1-status[i])
      ll[i] <- aux1*aux2
    }
  }
  #
  lll <- log(as.numeric(ll))
  return(sum(lll))
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
      Aux4 <- matrix(c(Aux3, Aux1, Aux2), (p+q+1),1)
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
      Aux4 <- matrix(c(Aux3, Aux1, Aux2), (p+q+1),1)
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
      Aux4 <- matrix(c(Aux3, Aux1, Aux2), (p+q+1),1)
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
      Aux4 <- matrix(c(Aux3, Aux1, Aux2), (p+q+1),1)
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
      Aux4 <- matrix(c(Aux3, Aux1, Aux2), (p+q+1),1)
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

  ep.aux <- sqrt(diag(solve(MI)))
  ep.beta <- ep.aux[1:p]
  ep.eta <- ep.aux[(p+1):(p+q)]

  if (dist != "exponential" && dist != "rayleigh") {
    ep.alpha <- ep.aux[(p+q+1)]
    ep.final <- c(ep.alpha, ep.beta, ep.eta)
  } else {
    ep.final <- c(ep.beta, ep.eta)
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
#
#
#
#' Fit a Mixture Rate Cure (MRC) Model with Custom Link Function
#'
#' This function fits a survival model to the provided data using a specified link function, with options for weights, iteration limits, and tolerance.
#'
#' @param time A numeric vector. The observed survival times.
#' @param status A numeric vector. The censoring indicator, where `1` typically indicates an event (e.g., death) and `0` indicates right censoring.
#' @param x A numeric matrix or data frame. Covariates for the survival model.
#' @param w A numeric vector. Weights for each observation. Defaults to equal weights if not provided.
#' @param maxit An integer. The maximum number of iterations for the optimization procedure. Defaults to `1000`.
#' @param tol A numeric value. The convergence tolerance for the optimization procedure. Defaults to `1e-5`.
#' @param link A character string. The link function to use. Defaults to `"logit"`. Options include `"logit"`, `"probit"`, `"plogit"`, `"rplogit"` or `"cauchit"`.
#' @param dist A character string. The distribution function to use. Defaults to `"weibull"`. Options include `"Exponential"`, `"Rayleigh"`, `"Log-normal"` or `"Log-logistic"`.
#' @param tau A numeric value. A threshold or cutoff parameter. Defaults to `1`.
#'
#' @return A list containing:
#'   \item{coefficients}{The estimated coefficients for the covariates.}
#'   \item{convergence}{A logical value indicating whether the optimization converged.}
#'   \item{iterations}{The number of iterations performed.}
#'   \item{logLik}{The log-likelihood of the fitted model.}
#'
#' @details
#' This function implements a custom survival model fitting procedure, allowing the user to specify a link function and other model parameters. It supports weighted observations and iterative optimization with user-defined tolerance and iteration limits.
#'
#' @examples
#' data(liver)
#' #head(liver)
#' model1 <- MCRfit(Surv(time,status)~age+medh+relapse+grade|-1+sex+age+medh+grade,data=liver)
#' summary(model1)
#' @export
MCRfit<-function(formula,data,dist="weibull",
                link="logit",tau=1,
                maxit = 1E3, tol = 1E-5){

  #

  x <- model.matrix(Formula(formula), data = data, rhs = 1)
  w <- model.matrix(Formula(formula), data = data, rhs = 2)

  mf <- model.frame(Formula(formula), data = data)
  model.aux <- model.response(mf)
  cc <- model.aux[, "status"]
  y <- model.aux[, "time"]

  # Initial estimates

  fit0 <- survreg(Surv(y, cc) ~ x - 1, dist = dist)
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
      aux<- exponential_sf(y,alpha=1,lambda) #para1[1]
    } else if(dist=="rayleigh"){
      aux<- rayleigh_sf(y,alpha=2,lambda) #para1[1]
    } else if(dist=="weibull"){
      aux<- weibull_sf(y,para1[1],lambda) #para1[1]
    } else if(dist=="lognormal"){
      aux<- lognormal_sf(y,para1[1],-log(lambda)) #para1[1]
    } else if(dist == "loglogistic"){
      aux <- loglogistic_sf(y, para1[1], lambda) #para1[1]
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

    fit <- survreg(Surv(y, cc) ~ x - 1,
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

    if (dist!="exponential" && dist!="rayleigh"){
      para2 <- c(alpha, beta, eta)
    }else{
      para2 <- c(beta, eta)
    }

    criteria <- (para2-para1)%*%(para2-para1)

    iter <- iter+1
    para1 <- para2
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

  n_alpha <- if (dist %in% c("weibull", "lognormal", "loglogistic")) 1 else 0
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

  if (!is.null(object$scale) && length(object$scale) == 1 && object$scale != 1) {
    std_alpha <- object$ep[1]  # assumindo que o erro padrão da escala é o primeiro
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

  cat("Coefficients (survival part):\n")
  print(format(x$coef_surv, digits = digits, nsmall = digits), quote = FALSE)

  cat("\nCoefficients (cure part):\n")
  print(format(x$coef_cure, digits = digits, nsmall = digits), quote = FALSE)

  if (!is.null(x$scale)) {
    cat("\nScale:\n")
    print(format(x$scale, digits = digits, nsmall = digits), quote = FALSE)
  }
}
