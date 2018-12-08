# load the source code of the functions to be tested

source("ARS.R")



# context with one test that groups expectations

context("Tests for ARS() and auxiliary functions")



# Tests for the overall function

# Normal distribution test

test_that("ars() used for normal dist", {

  g <- function(x){return(dnorm(x,0,1))}

  sample_norm <- ars(g,1000,-3,3) # Generate normal samples from ars()

  sample_norm_cdf <- c()

  true_norm_cdf <- c()

  test_norm_interval <- seq(-3,3,0.5)# Test array

  for (i in 1:(length(test_norm_interval)-1)){

    # The difference of sample normal interval cdf

    sample_norm_cdf[i] <- (sum(sample_norm < test_norm_interval[i+1])-sum(sample_norm < test_norm_interval[i]))/1000

    # The difference of true normal interval cdf 

    true_norm_cdf[i] <- (pnorm(test_norm_interval[i+1])-pnorm(test_norm_interval[i]))/(pnorm(3)-pnorm(-3))

  }

  # Take the maximum error

  max_norm_error <- max(abs(sample_norm_cdf-true_norm_cdf))

  

  expect_lt(max_norm_error, 0.1)# Max_error should be less than the threshold

})



# Uniform distribution test

test_that("ars() used for uniform dist", {

  g <- function(x){return(dunif(x, min = 1, max = 10))}

  sample_unif <- ars(g,1000,2,9) # Generate uniform samples from ars()

  sample_unif_cdf <- c()

  true_unif_cdf <- c()

  test_unif_interval <- seq(2,9,0.5)# Test array

  for (i in 1:(length(test_unif_interval)-1)){

    # The difference of sample uniform interval cdf

    sample_unif_cdf[i] <- (sum(sample_unif < test_unif_interval[i+1])-sum(sample_unif < test_unif_interval[i]))/1000

    # The difference of true uniform interval cdf 

    true_unif_cdf[i] <- (punif(test_unif_interval[i+1],1,10)-punif(test_unif_interval[i],1,10))/(punif(9,1,10)-punif(2,1,10))

  }

  # Take the maximum error

  max_unif_error <- max(abs(sample_unif_cdf-true_unif_cdf))

  

  expect_lt(max_unif_error, 0.01)# Max_error should be less than the threshold

})





# Exponential distribution test

#sample_exp <- ars("dexp(x, rate = 1, log = FALSE)",1000,1,10)



# Gamma distribution test

test_that("ars() used for gamma dist", {

  g <- function(x){return(dgamma(x, 2, rate = 1))}

  sample_gamma <- ars(g,1000,1,5) # Generate gamma samples from ars()

  sample_gamma_cdf <- c()

  true_gamma_cdf <- c()

  test_gamma_interval <- seq(1,10,0.5)# Test array

  for (i in 1:(length(test_gamma_interval)-1)){

    # The difference of sample gamma interval cdf

    sample_gamma_cdf[i] <- (sum(sample_gamma < test_gamma_interval[i+1])-sum(sample_gamma < test_gamma_interval[i]))/1000

    # The difference of true gamma interval cdf 

    true_gamma_cdf[i] <- (pgamma(test_gamma_interval[i+1],2,1)-pgamma(test_gamma_interval[i],2,1))/(pgamma(5,2,1)-pgamma(1,2,1))

  }

  # Take the maximum error

  max_gamma_error <- max(abs(sample_gamma_cdf-true_gamma_cdf))

  

  expect_lt(max_gamma_error, 0.01)# Max_error should be less than the threshold

})





# Unit tests



# Intersect_Z(X_k, h, lb, ub)

test_that("Intersect_Z fuction", {

  X_k <- c(-3,0,3)

  h <- function(x){return(log(dnorm(x)))}

  Z_k <- Intersect_Z(c(-3,0,3), h, lb=-3, ub=3)

  

  expect_lt(Z_k[2], X_k[2])# Z_k[2] < X_k[2]

  expect_lt(Z_k[3], X_k[3])# Z_k[3] < X_k[3]

  expect_gt(Z_k[2], X_k[1])# Z_k[2] > X_k[1]

  expect_length(Z_k, 4)# There should be 4 numbers in Z_k for this example

 })



# U_k(x, h, X_k, Z_k)

test_that("U_k fuction", {

  X_k <- c(-3,0,3)

  h <- function(x){return(log(dnorm(x)))}

  Z_k <- Intersect_Z(c(-3,0,3), h, lb=-3, ub=3)

  # Two methods to get j, j_1 is used in the package

  j_1 <- min(which(0 < Z_k))-1

  j_2 <- max(which(0 > Z_k))

  

  expect_equal(j_1, j_2)# j_1 = j_2: Use another methhod to find j(j_2), should be equal to j_1

  # j_1 should be a double number

  expect_type(j_1, 'double')

})



# Integrate_ExpU_k(h, X_k, Z_k)

test_that("Integrate_ExpU_k", {

  X_k <- c(-3,0,3)

  h <- function(x){return(log(dnorm(x)))}

  Z_k <- Intersect_Z(c(-3,0,3), h, lb=-3, ub=3)

  # Two methods to get integral from Z_k[1] to Z_k[2]

  # exp_u_k_1 and integral_1 are used in the package because of for loop

  exp_u_k_1 <- function(x,j){

    return(exp(h(X_k[j]) + (x - X_k[j]) * gradient(h, X_k[j])[1][1]))

  }

  integral_1 <- integrate(f=exp_u_k_1,lower=Z_k[1], upper=Z_k[2], j=1)$value

  exp_u_k_2 <- function(x){

    return(exp(h(X_k[1]) + (x - X_k[1]) * gradient(h, X_k[1])[1][1]))

  }

  integral_2 <- integrate(f=exp_u_k_2,lower=Z_k[1], upper=Z_k[2])$value

  expect_equal(integral_1, integral_2)# integral_1 = integral_2

  # integral_1 should be of type double

  expect_type(integral_1, 'double')

})



# L_k(x, h, X_k)

test_that("L_k function", {

  X_k <- c(-3,0,3)

  h <- function(x){return(log(dnorm(x)))}

  Z_k <- Intersect_Z(c(-3,0,3), h, lb=-3, ub=3)

  # Two methods to get j, j_1 is used in the package

  j_1 <- max(which(0.5 > X_k))

  j_2 <- min(which(0.5 < X_k))-1

  # 5 is not in the range of X_k in this example, L_k should be -Inf

  L <- L_k(5,h,X_k)

  expect_equal(j_1, j_2)# j_1 = j_2: Use another methhod to find j(j_2), 

  #should be equal to j_1

  # j_1 should be a single integer

  expect_type(j_1, 'integer')

  expect_equal(exp(L), 0)

})



# Module tests



# Initial_X(h,lb,ub)

test_that("Initial_X step module", {

  X_k <- c(-3,0,3)

  h <- function(x){return(log(dnorm(x)))}

  Z_k <- Intersect_Z(c(-3,0,3), h, lb=-3, ub=3)

  # Get the maximizer for h(x) function

  xmax <- optimize(f = h, interval = c(-3, 1), lower = -3, upper = 1, maximum = TRUE)$maximum

  if ((abs(xmax - -3) < 1e-8) | (abs(xmax - 3) < 1e-8)) {

    Xinit <- c(-3,(-3+1)/2,3)

  }

  else {

    Xinit <- c(-3,xmax,3)

  }

  expect_equal(xmax, 0)# In the example for normal distribution, xmax should be 0.

  # Xinit[2] should be 0 not -1

  expect_equal(Xinit[2], 0)

})



# Sample_step(u, S_k, lb, ub)



test_that("Sample_step Module", {

  # Initialization

  lb <- -3

  ub <- 3

  # Suppose s(x) as normal dist pdf

  g <- function(x){return(dnorm(x,0,1))}

  # Calculate cdf at the boundary of x

  cdf <- function(x) integrate(g, lb, x)$value

  # 0.9772499 is pnorm(2,0,1)

  f <- function(x) ((cdf(x)/cdf(ub)) -  0.9772499)

  xstar <- uniroot(f, lower=lb, upper=ub)$root

  

  expect_lt((xstar - 2), 1e-5)# xstar should be really close to 2

})



# Rejection_step(x_star,h,X_k,Z_k)

test_that("Rejection_step module", {

  X_k <- c(-3,0,3)

  h <- function(x){return(log(dnorm(x)))}

  Z_k <- Intersect_Z(c(-3,0,3), h, lb=-3, ub=3)

  set.seed(1)

  w = runif(1)

  # In the setup, w=0.27, x_star=2(in the range of X_k)

  # So threshold_1=0.22, not satisfied

  # threshold_2=0.606, satisfied

  

  threshold_1 = exp(L_k(2, h, X_k) - U_k(2, h, X_k, Z_k))

  threshold_2 = exp(h(2) - U_k(2, h, X_k, Z_k))

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

  accept_update <- c(accept, update)

  

  # accept_update should be (TRUE,TRUE) in this example

  expect_equal(accept_update[1], TRUE)

  expect_equal(accept_update[2], TRUE)

})
