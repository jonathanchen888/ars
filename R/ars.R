library(assertthat)
library(rootSolve)
## Equation 1: Compute intersection points Z
Intersect_Z <- function(X_k, h, lb, ub){

  # X_k is a list with points,
  # h is a concave and differentiable function
  # lb is the lower bound of domain X
  # ub is the upper bound of domain X

  Z <- rep(0, length(X_k)-1)

  ## loop through 1:k-1 for Z
  for (j in 1:(length(X_k)-1)){
    x_j1 <- X_k[j+1]
    x_j <- X_k[j]
    Numerator <- h(x_j1) - h(x_j) - x_j1*gradient(h,x_j1)[1][1] + x_j*gradient(h,x_j)[1][1]
    Denominator <- gradient(h,x_j)[1][1] - gradient(h,x_j1)[1][1]
    Z[j] <- Numerator/(Denominator+1e-8)
  }

  Z <- c(lb, Z, ub)

  return(Z)
}


## Equation 2: Upper bound, U_k(x)
U_k <- function(x, h, X_k, Z_k){

  # x is a point
  # h is a concave and differentiable function
  # X_k is a list with points,
  # Z_k is a list get from Intersect_Z

  j <- min(which(x < Z_k))-1  ## j in [1:k]
  return(h(X_k[j]) + (x - X_k[j]) * gradient(h, X_k[j])[1][1])
}
# U_k(0, h, X_k, Z_k)


Integrate_ExpU_k <- function(h, X_k, Z_k){
  integ_sum <- c()
  exp_u_k <- function(x,j){
    return(exp(h(X_k[j]) + (x - X_k[j]) * gradient(h, X_k[j])[1][1]))
  }
  for (i in (2:length(Z_k))){
    integ <- integrate(f=exp_u_k, lower=Z_k[i-1], upper=Z_k[i], j=i-1)$value
    integ_sum = c(integ_sum,integ)
  }
  return(sum(unlist(integ_sum)))
}
# integ_sum <- Integrate_ExpU_k(h, X_k, Z_k)


## Equation 3: Envelope_Function
Envelope_Function <- function(x, h, X_k, Z_k, lb, ub, integ_sum){
  exp_Uk <- exp(U_k(x, h, X_k, Z_k))
  return(exp_Uk/(integ_sum+1e-8))
}

## Equation 4: Lower bound, L_k(x)
L_k <- function(x, h, X_k){

  # x is a point
  # h is a concave and differentiable function
  # X_k is a list with points

  if(x < min(X_k) | x > max(X_k)){
    return(-Inf)
  }else{
    j <- max(which(x >= X_k))
    Numerator <- (X_k[j+1]-x)*h(X_k[j]) + (x-X_k[j])*h(X_k[j+1])
    Denominator <- X_k[j+1] - X_k[j]
    return(Numerator/(Denominator+1e-8))
  }
}


## Generate a sample from envelope S_k(x) of section 2.2.2
Sample_step <- function(u, S_k, lb, ub){
  cdf <- function(x) integrate(S_k, lb, x)$value
  f <- function(x) ((cdf(x)/cdf(ub)) - u)
  x_star <- uniroot(f, lower=lb, upper=ub)$root
  # cdf <- function(x) integrate(S_k, lb, x)$value
  # qdf <- function(x) optimize(function(z)(cdf(z)-x)^2,c(lb,ub))$minimum
  # rdf<-function(n) sapply(runif(n),qdf)
  # x_star <- qdf(u)

  return(x_star)
}

## Accept/Rejection test of section 2.2.2
Rejection_step <- function(x_star,h,X_k,Z_k){
  w = runif(1)

  # Rejection tests
  threshold_1 = exp(L_k(x_star, h, X_k) - U_k(x_star, h, X_k, Z_k))
  threshold_2 = exp(h(x_star) - U_k(x_star, h, X_k, Z_k))
  if( w <= threshold_1){
    accept = TRUE
    update = FALSE
  }
  else if(w <= threshold_2){
    update = TRUE
    accept = TRUE
  }
  else{
    update = TRUE
    accept = FALSE
  }

  ## whether to accept candidate sample point, and whether add the point into X
  return(c(acc=accept,upd=update))
}


Initial_X <- function(h,lb,ub) {

  # h is a concave and differentiable function
  # lb is the lower bound of domain X
  # ub is the upper bound of domain X

  x_max <- optimize(f = h, interval = c(lb, ub), lower = lb, upper = ub, maximum = TRUE)$maximum

  ## if max is at upper bound or lower bound
  if ((abs(x_max - ub) < 1e-8) | (abs(x_max - lb) < 1e-8)) {
    X_init <- c(lb,(lb+ub)/2,ub)
  }
  else {
    X_init <- c(lb,x_max,ub)
  }

  return(X_init)
}

#' Adaptive Rejection Sampling
#'
#' Based on the Gilks & Wild (1992) paper
#'
#' @param g a pdf (probability density function)
#' @param n number of samples to generate
#' @param lb lower bound of domain
#' @param ub upper bound of domain
#' @return samples of length n from g
#' @export

##########################################################################
ars <- function(g, n, lb, ub){
  ## Some simple inputs tests
  assert_that(is.numeric(n))
  assert_that(is.numeric(lb))
  assert_that(is.numeric(ub))
  if (lb >= ub)
    stop("Not valid domain.")
  h <- function(x){return(log(g(x)))}

  ## Test for log concave, lower bound & upper bound derivative
  if (!(gradient(h,lb)[1][1]>=0))
    stop("The gradient of the log density at lower bound should be bigger than 0")
  if (!(gradient(h,ub)[1][1]<=0))
    stop("The gradient of the log density at upper bound should be less than 0")
  test_array <- seq(lb,ub,0.01)
  for (i in 1:(length(test_array)-1)){
    gradient_diff <- gradient(h,test_array[i])[1][1] - gradient(h,test_array[i+1])[1][1]
    if(gradient_diff < 0)
      stop("The density function is not log concave")
  }


  ## Initial X_k, batchsize and samples
  X_init <- Initial_X(h,lb,ub)
  X_k <- X_init
  samples = c()
  batchsize <- 1.2*n  ## batchsize to sample

  while(length(samples) < n){

    # Get intersection points
    Z_k <- Intersect_Z(X_k,h,lb,ub)
    integ_sum <- Integrate_ExpU_k(h, X_k, Z_k)

    ## Sampling step
    X_unif <- runif(batchsize)
    S_k <- function(x){return(Envelope_Function(x, h, X_k, Z_k, lb, ub, integ_sum))}
    X_stars <- sapply(X_unif, Sample_step, S_k=S_k, lb=lb, ub=ub)
    out <- sapply(X_stars, Rejection_step, h=h, X_k=X_k, Z_k=Z_k)
    X_accepts <- out[1,]
    X_to_adds <- out[2,]
    ## Updating step
    samples = c(samples, X_stars[X_accepts=TRUE])
    X_k = sort(c(X_k, X_stars[X_to_adds=TRUE]))

    ## Update batchsize for next loop
    batchsize <- n-length(samples)
  }

  return(samples[1:n])
}
