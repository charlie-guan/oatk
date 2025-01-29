#' Ridge regression
#'
#' @description
#' Function for conducting ridge and for finding lambda using fast leave-one-out cross-validation (LOOCV) if
#' lambda is a vector.
#'
#' @param Y Response variable
#' @param U U in the SVD X = UDV^t,
#' @param D D in the SVD X = UDV^t,
#' @param V V in the SVD X = UDV^t,
#' @param lam Ridge parameter. If lam is a scalar, ridge will be conducted with that lam.
# If lam is a vector, LOOCV is used to find the best value. WARNING: In degenerate cases
#' (X nearly singular & training fit nearly perfect),
# neither LOOCV nor GCV are reliable and select lambda way too small. For any such
# examples, you should use K-fold CV with K ~ 10 or so to select lambda.
#'
#' @return A list containing beta_hat (the regression coefficients) and lam.ridge (the corresponding ridge parameter)
ridge.reg <- function(Y, U, D, V, lam) {
  m = length(lam)
  n = length(Y)
  if (m > 1) {
    MSE.CV <- rep(0,m)
    #    MSE.GCV <- rep(0,m) #In case we want to use GCV instead of LOOCV
    for (i in 1:m) {
      d.sq.lam <- D^2/(D^2+lam[i])
      Y.hat <- U%*%(d.sq.lam*(t(U)%*%Y)) #ridge regression fitted response values
      e <- as.numeric(Y-Y.hat) #ridge regression residual errors
      PRESS.den <- apply(U, 1, function(x) 1-sum(x^2*d.sq.lam)) #1-H_{i,i} for PRESS denominator
      MSE.CV[i] <- sum((e/PRESS.den)^2)/n
      #      GCV.den <- 1-sum(apply(U, 1, function(x) sum(x^2*d.sq.lam)))/n #1-trace(H)/n for GCV denominator
      #      MSE.GCV[i] <- sum(e^2)/(n*GCV.den)
    } #end of for loop
    #    plot(lam,MSE.CV)
    lam.best <-lam[which.min(MSE.CV)]
  }else {
    lam.best <-lam
  } #end of if-else
  beta.hat <- V%*%((D/(D^2+lam.best))*(t(U)%*%Y)) #ridge regression parameters
  return( list(beta.hat = beta.hat, lam.best = lam.best))
} #end of ridge.lam function


#' L2 norm
#'
#' @description
#' Calculates L2 norm of a vector
#'
#' @param x numeric vector
#'
#' @return A scalar corresponding to the L2 norm
l2_norm = function(x) {
  return(sqrt(sum(x^2)))
}


#' eBH
#'
#' @description This conducts the eBH procedure, proposed by Wang and Ramdas (2022),
#' to conduct variable selection with controlled false discovery rate (FDR) using e-values.
#'
#' @param e_values Numeric vector containing e-values
#' @param alpha The desired FDR level.
#'
#' @return A numeric vector corresponding to the rejection set.
eBH = function(e_values, alpha) {
  o = order(e_values, decreasing=TRUE)
  p = length(e_values)

  bh = max(which(e_values[o] * (1:p) / p >= 1 / alpha), 0)
  if (bh==0) {
    rej = NULL
  } else {
    rej = o[1:bh]
  }

  return(rej)
}
