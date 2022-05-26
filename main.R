####################################################
# Here is the R code for the maximum likelihood estimation 
# of alpha-stable Double Autoregressive Models.

#' If you use sDAR model in your research, please cite our work
#' by using the following BibTeX entry:
#' @article{sDAR2022,
#'   title={Maximum likelihood estimation for $\alpha$-stable double autoregressive models},
#'   author={Li, Dong and Tao, Yuxin and Yang, Yaxin and Zhang, Rongmao},
#'   journal={submitted},
#'   year={2022+}
#' }

####################################################

library(stabledist)
library(MASS)
source('stable_density.R')
source('stable_dx.R')
source('stable_dalpha.R')
load("Ef_bA_estimate.RData")

####################################################
# Calculate the Lyapunov exponent gamma
####################################################
Lyapunov <- function(a,b,alpha){
  gamma <- .integrate2(function(x) log(abs(a^2-b*x^2))*dstable(x,alpha=alpha,beta=0),
                       0,Inf,subdivisions = 2000,rel.tol=.Machine$double.eps^.05,abs.tol=.Machine$double.eps^.05)
  return(gamma)
}

####################################################
# Maximum likelihood estimation of sDAR(1) model
####################################################

### Log-likelihood function
lik <- function(theta){
  a <- theta[1]; b <- theta[2]; c <- theta[3]; alpha <- theta[4]
  sum <- 0
  for (t in 2:(n+1)){
    f <- dstable((y[t]-a*y[t-1])/sqrt(c+b*(y[t-1]^2)), alpha=alpha, beta=0)
    sum <- sum-log(c+b*(y[t-1]^2))/2+log(f)
  }
  return(-sum)
}
### Gradient function 
grr <- function(theta){
  a <- theta[1]; b <- theta[2]; c <- theta[3]; alpha <- theta[4]
  sum <- rep(0,4)
  for (t in 2:(n+1)){
    qt <- (y[t]-a*y[t-1])/sqrt(c+b*(y[t-1]^2))
    f1 <- dstable_dx_new(qt,alpha)
    f2 <- dstable(qt, alpha=alpha, beta=0)
    deriv <- f1 / f2
    da <- -deriv*y[t-1]/sqrt(c+b*(y[t-1]^2))
    dc <- -1/2/(c+b*(y[t-1]^2))*(1+deriv*qt)
    db <- dc*(y[t-1]^2)
    f3 <- dstable_da_new(qt,alpha)
    dalpha <- f3 / f2
    sum <- sum + c(da,db,dc,dalpha)
  }
  return(-sum)
}

#### Log-likelihood function for diagnostic checking
lik1 <- function(theta){
  a <- theta[1]; b <- theta[2]; c <- theta[3]
  sum <- 0
  for (t in 2:(n+1)){
    f <- dstable((y[t]-a*y[t-1])/sqrt(c+b*(y[t-1]^2)), alpha=alpha, beta=0)
    sum <- sum-log(c+b*(y[t-1]^2))/2+log(f)
  }
  if (-sum==Inf) {sum <- -1e10}
  return(-sum)
}
### Gradient function for diagnostic checking
grr1 <- function(theta){
  a <- theta[1]; b <- theta[2]; c <- theta[3]
  sum <- rep(0,3)
  for (t in 2:(n+1)){
    qt <- (y[t]-a*y[t-1])/sqrt(c+b*(y[t-1]^2))
    f1 <- dstable_dx_new(qt,alpha)
    f2 <- dstable(qt,alpha=alpha,beta=0)
    deriv <- f1 / f2
    da <- -deriv*y[t-1]/sqrt(c+b*(y[t-1]^2))
    dc <- -1/2/(c+b*(y[t-1]^2))*(1+deriv*qt)
    db <- dc*(y[t-1]^2)
    sum <- sum + c(da,db,dc)
  }
  return(-sum)
}


### Maximum Likelihood Estimation and Diagnostic Checking
sDAR <- function(y, a_0=0.5, b_0=0.5, c_0=0.5, alpha_0=1){
  ### Use the first value as the initial value y_0
  n <- length(y)-1
  theta <- try(constrOptim(c(a_0,b_0,c_0,alpha_0), f=lik, grad=grr,
                         ui = rbind(c(0,1,0,0),c(0,0,1,0),c(0,0,0,1),c(0,0,0,-1)),
                         ci = c(0,0,0,-2))$par
  )  
  a <- theta[1]; b <- theta[2]; c <- theta[3]; alpha <<- theta[4]
  
  ### Calculating the estimated standard deviations and Lyapunov exponents
  gamma_1 <- Lyapunov(a,b,alpha)
  ASD <- rep(0,16); SD_1 <- rep(NA,4); SD_2 <- rep(NA,4); gamma_2 <- 0
  if (gamma_1 < 0){### Stationary case
    ASD[1] <- y[1]^2/(c+b*y[1]^2)
    ASD[2] <- (y[1]^2/(c+b*y[1]^2))^2
    ASD[3] <- y[1]^2/(c+b*y[1]^2)^2
    ASD[4] <- y[1]^2/(c+b*y[1]^2)
    ASD[5] <- 1/(c+b*y[1]^2)^2
    ASD[6] <- 1/(c+b*y[1]^2)
    for (t in 2:(n+1)){
      eta_t <- (y[t] - a * y[t-1])/sqrt(c + b * (y[t-1]^2))
      ASD[1] <- ASD[1] + y[t]^2/(c+b*y[t]^2)
      ASD[2] <- ASD[2] + (y[t]^2/(c+b*y[t]^2))^2
      ASD[3] <- ASD[3] + y[t]^2/(c+b*y[t]^2)^2
      ASD[4] <- ASD[4] + y[t]^2/(c+b*y[t]^2)
      ASD[5] <- ASD[5] + 1/(c+b*y[t]^2)^2
      ASD[6] <- ASD[6] + 1/(c+b*y[t]^2)
      deriv_log_x <- dstable_dx_new(eta_t,alpha)/dstable(eta_t, alpha=alpha, beta=0)
      deriv_log_A <- dstable_da_new(eta_t,alpha)/dstable(eta_t, alpha=alpha, beta=0)
      ASD[7] <- ASD[7] + deriv_log_x^2
      ASD[8] <- ASD[8] + (deriv_log_x^2) * (eta_t^2)
      ASD[9] <- ASD[9] + deriv_log_A^2
      ASD[10] <- ASD[10] + deriv_log_x * deriv_log_A * eta_t
      gamma_2 <- gamma_2 + log(abs(a+eta_t*sqrt(b)))
    }
    ASD[c(1:6)] <- ASD[c(1:6)]/(n+1)
    ASD[c(7:10)] <- ASD[c(7:10)]/n
    gamma_2 <- gamma_2/n
    
    ASD[11] <- ASD[1]*ASD[7]
    ASD[12] <- 0.25*ASD[2]*(ASD[8]-1)
    ASD[13] <- 0.25*ASD[3]*(ASD[8]-1)
    ASD[14] <- -0.5*ASD[4]*ASD[10]
    ASD[15] <- 0.25*ASD[5]*(ASD[8]-1)
    ASD[16] <- -0.5*ASD[6]*ASD[10]
    sigma <- matrix(c(ASD[11],0,0,0,
                     0,ASD[12],ASD[13],ASD[14],
                     0,ASD[13],ASD[15],ASD[16],
                     0,ASD[14],ASD[16],ASD[9]),4,4)
    d <- try(solve(sigma))
    if (!'try-error' %in% class(d)) {SD_1 <- sqrt(diag(d/n))}
    
    Ef_aa <- .integrate2(function(x) 2*dstable_dx_new(x,alpha)^2/dstable(x,alpha=alpha,beta=0),
                           0,Inf,subdivisions = 2000,rel.tol =.Machine$double.eps^.05,abs.tol =.Machine$double.eps^.05)
    Ef_bb <- .integrate2(function(x) 2*x^2*dstable_dx_new(x,alpha)^2/dstable(x,alpha=alpha,beta=0),
                           0,Inf,subdivisions = 2000,rel.tol=.Machine$double.eps^.05,abs.tol =.Machine$double.eps^.05)
    Ef_alpha <- .integrate2(function(eta) 2*dstable_da_new(eta,alpha)^2/dstable(eta,alpha=alpha,beta=0),
                              0,Inf,subdivisions = 2000,rel.tol=.Machine$double.eps^.05,abs.tol =.Machine$double.eps^.05)
    if (alpha<0.76) { ### numerical issues ---> replace it by estimated values already calculated in result2
      Ef_bA <- result2[2,round(alpha*100)-29]
    } else {
      Ef_bA <- .integrate2(function(x) 2*x*dstable_dx_new(x,alpha)/dstable(x,alpha=alpha,beta=0)*dstable_da_new(x,alpha),
                             0,Inf,subdivisions = 2000,rel.tol=.Machine$double.eps^.05,abs.tol =.Machine$double.eps^.05)
    }
    sig_aa <- ASD[1]*Ef_aa
    sig_bb <- 0.25*ASD[2]*(Ef_bb-1)
    sig_bc <- 0.25*ASD[3]*(Ef_bb-1)
    sig_bA <- -0.5*ASD[4]*Ef_bA
    sig_cc <- 0.25*ASD[5]*(Ef_bb-1)
    sig_cA <- -0.5*ASD[6]*Ef_bA
    sig_AA <- Ef_alpha
    sigma <- matrix(c(sig_aa,0,0,0,
                      0,sig_bb,sig_bc,sig_bA,
                      0,sig_bc,sig_cc,sig_cA,
                      0,sig_bA,sig_cA,sig_AA),4,4)
    d <- try(solve(sigma))
    if (!'try-error' %in% class(d)) {SD_2 <- sqrt(diag(d/n))}
  } else {### Explosive case
    for (t in 2:(n+1)){
      eta_t <- (y[t] - a * y[t-1])/sqrt(c + b * (y[t-1]^2))
      deriv_log_x <- dstable_dx_new(eta_t,alpha)/dstable(eta_t, alpha=alpha, beta=0)
      deriv_log_A <- dstable_da_new(eta_t,alpha)/dstable(eta_t, alpha=alpha, beta=0)
      ASD[7] <- ASD[7] + deriv_log_x^2
      ASD[8] <- ASD[8] + (deriv_log_x^2) * (eta_t^2)
      ASD[9] <- ASD[9] + deriv_log_A^2
      ASD[10] <- ASD[10] + deriv_log_x * deriv_log_A * eta_t
      gamma_2 <- gamma_2 + log(abs(a + eta_t * sqrt(b)))
    }
    ASD[c(7:10)] <- ASD[c(7:10)]/n
    gamma_2 <- gamma_2/n
    
    ASD[11] <- 1/b * ASD[7]
    ASD[12] <- 1/(4 * b^2) * (ASD[8]-1)
    ASD[14] <- -1/(2 * b) * ASD[10]
    sigma <- matrix(c(ASD[11],0,0,
                      0,ASD[12],ASD[14],
                      0,ASD[14],ASD[9]),3,3)
    d <- try(solve(sigma))
    if (!'try-error' %in% class(d)) {SD_1[c(1,2,4)] <- sqrt(diag(d/n))}
    
    Ef_aa <- .integrate2(function(x) 2*dstable_dx_new(x,alpha)^2/dstable(x,alpha=alpha,beta=0),
                           0,Inf,subdivisions = 2000,rel.tol =.Machine$double.eps^.05,abs.tol =.Machine$double.eps^.05)
    Ef_bb <- .integrate2(function(x) 2*x^2*dstable_dx_new(x,alpha)^2/dstable(x,alpha=alpha,beta=0),
                           0,Inf, subdivisions = 2000,rel.tol=.Machine$double.eps^.05,abs.tol =.Machine$double.eps^.05)
    Ef_alpha <- .integrate2(function(eta) 2*dstable_da_new(eta,alpha)^2/dstable(eta,alpha=alpha,beta=0),
                              0,Inf,subdivisions = 2000,rel.tol=.Machine$double.eps^.05,abs.tol =.Machine$double.eps^.05)
    if (alpha<0.76) {
      Ef_bA <- result2[2,round(alpha*100)-29]
    } else {
      Ef_bA <- .integrate2(function(x) 2*x*dstable_dx_new(x,alpha)/dstable(x,alpha=alpha,beta=0)*dstable_da_new(x,alpha),
                           0,Inf,subdivisions = 2000,rel.tol=.Machine$double.eps^.05,abs.tol =.Machine$double.eps^.05)
    }
    sig_aa <- 1/b * Ef_aa
    sig_bb <- 1/(4*b^2) * (Ef_bb-1)
    sig_bA <- -1/(2*b) * Ef_bA
    sig_AA <- Ef_alpha
    sigma <- matrix(c(sig_aa,0,0,
                      0,sig_bb,sig_bA,
                      0,sig_bA,sig_AA),3,3)
    d <- try(solve(sigma))
    if (!'try-error' %in% class(d)) {SD_2[c(1,2,4)] <- sqrt(diag(d/n))}
  }
  
  ### Diagnostic checking for the residuals
  esti <- try(constrOptim(c(a,b,c), f=lik1, grad=grr1,
                          ui = rbind(c(0,1,0),c(0,0,1)),
                          ci = c(0,0))$par
  )
  a <- esti[1]; b <- esti[2]; c <- esti[3]
  res <- (y[2:(n+1)] - a * y[1:n])/sqrt(c + b * (y[1:n]^2))
  
  U <- pstable(res, alpha=alpha, beta=0)
  U_sort <- sort(U); res_sort <- sort(res)
  nu <- c(U_sort,1); m <- length(res)
  g_dot <- matrix(0,m,2)
  for (j in 1:m){
    g_dot[j,1] <- 1
    F_inverse <- res_sort[j]
    g_dot[j,2] <- dstable_dx_new(F_inverse,alpha=alpha)/dstable(F_inverse,alpha=alpha,beta=0)
  }
  Ck <- matrix(0,m,4)
  Ck[m,] <- c(g_dot[m,] %*% t(g_dot[m,]) * (nu[m+1] - nu[m]))
  Dk <- matrix(0,m,2)
  Dk[m,] <- g_dot[m,]
  for (i in (m-1):1){
    Ck[i,] <- Ck[i+1,] + c(g_dot[i,] %*% t(g_dot[i,]) * (nu[i+1] - nu[i]))
    Dk[i,] <- Dk[i+1,] + g_dot[i,]
  }
  W <- NULL; add <- t(g_dot[1,]) %*% solve(matrix(Ck[1,],2,2)) %*% Dk[1,] * nu[1]
  W[1] <- 1/m - add/m
  for (j in 2:m){
    add <- try(add + t(g_dot[j,]) %*% solve(matrix(Ck[j,],2,2)) %*% Dk[j,] * (nu[j]-nu[j-1]))
    if ('try-error' %in% class(add)) {
      break
    } else {
      W[j] <- j/m-add/m
    }
  }
  Wn_hat <- max(abs(W))*sqrt(m)
  diagnos <- Wn_hat > c(1.9600, 2.2414, 2.8070)
  diagnos_result <- c("Accept the null", "Reject the null")
  
  ### Summarize the result
  MLE <- matrix(0,3,4, dimnames = list(c("Estimation","ASD_int","ASD_res"),
                                       c("a","b","c","alpha")))
  MLE[1,] <- theta; MLE[2,] <- SD_1; MLE[3,] <- SD_2
  result <- list("Maximum likelihood estimation" = round(MLE,4))
  result$"Lyapunov_int (recommended)" <- gamma_1
  result$"Lyapunov_res" <- gamma_2

  ### Goodness-of-fit
  good <- matrix(0,4,1, dimnames = list(c("log-likelihood","MSE","AIC","BIC"),
                                        c("Goodness-of-fit")))
  good[1,1] <- -lik(theta)
  good[2,1] <- mean((y[2:n]-a*y[1:(n-1)])^2)
  good[3,1] <- 2*4-2*(-lik(theta))
  good[4,1] <- 2*log(n)-2*(-lik(theta))
  result$"Goodness-of-fit" <- good
  result$"Diagnostic checking with null hypothesis: alpha_0=alpha_hat" <- 
    matrix(diagnos_result[diagnos+1],3,1,dimnames = 
             list(c("significance level of 10%: ","significance level of 5%: ",
                    "significance level of 1%: "), paste("Null: alpha_0 =",round(alpha,4))))
  return(result)
}


####################################################
# demo: a toy example
####################################################
### Simulate a sDAR series y with true value theta0 = c(a0,b0,c0,alpha0)
### length n+1 with initial value 0
a0 <- 1; b0 <- 0.3; c0 <- 0.5; alpha0 <- 1.5
n <- 400; flag <- TRUE
while (flag){ ### In case the series go beyond the numerical upper bound
  y <- c(0)
  for (t in 2:(n+1)){
    eta_t <- rstable(1, alpha = alpha0, beta=0)
    y[t] <- a0 * y[t-1] + eta_t * sqrt(c0 + b0 * (y[t-1]^2))
  }
  if (max(abs(y[!is.na(y)]))!=Inf) {flag <- FALSE}
}
summary(y); plot(y, type="l")

### Initial guess of parameters (for example)
ly <- y[2:n]/sqrt(1+y[1:(n-1)]^2); lx <- y[1:(n-1)]/sqrt(1+y[1:(n-1)]^2)
lmod <- rlm(ly ~ lx); a_0 <- as.numeric(lmod$coefficients[2])
b_0 <- c_0 <- 0.5; alpha_0 <- 1

### sDAR(1) model fitting
sDAR_esti <- sDAR(y,a_0,b_0,c_0,alpha_0)
sDAR_esti


