library(SNPknock)


#' Benjamini-Hochberg procedure
#'
#' @description
#' This function conducts the Benjamini-Hochberg procedure, which controls
#' the false discovery rate (FDR) in multiple hypothesis tests. The procedure is based on Benjamini and Hochberg (1995)
#'
#' @param y Response variable
#' @param X Design matrix
#' @param q The desired FDR level
#'
#' @return A list containing the vector rej, which are the indices of the columns of X corresponding to the non-null variables.
#'
#' @importFrom stats lm
#' @export
bh = function(y, X, q=0.1) {
  p = ncol(X)

  # obtain p-values of regression coefficients
  mdl = lm(Y ~ .-1, data = data.frame("Y" = y, X))
    p_values = summary(mdl)$coefficients[, 4]

  # run the BH procedure
  o = order(p_values)
  bh = max(which(p_values[o] <= 1:p * q / p), 0)
  if (bh==0) rej = integer() else rej = o[1:bh]

  return(list(rej=rej))
}



#' The Knockoff Filter
#'
#' @description This is a variable selection procedure that controls the false discovery rate, as proposed by Barber and Candes (2015)
#'
#' @param y Response variable
#' @param X Design matrix
#' @param q The desired FDR level
#' @param offset Offset parameter when estimated the FDR. Setting c=1 results in a more conservative FDR result.
#'
#' @return A list containing the vector rej, which are the indices of the columns of X corresponding to the non-null variables,
#' the vector W which are the test statistics for each variable, and the numeric threshold used to create the rejection set.
#'
#' @importFrom knockoff create.fixed stat.lasso_lambdasmax knockoff.threshold
#'
#' @export
bc = function(y, X, q=0.1, offset=0) {
  Xk <- create.fixed(X)$Xk
  W <- stat.lasso_lambdasmax(X, Xk, y)
  threshold = knockoff.threshold(W, fdr=q, offset=offset)
  rej = which(W >= threshold)

  return(list(W=W, rej=rej, threshold=threshold))
}



#' knockoffDMC
#'
#' @description This is a knockoff procedure specifically made for design matrix that is generated
#' from discrete Markov chains, developed by Sesia et al. (2019)
#'
#' @param y Response variable
#' @param X Design matrix
#' @param q The desired FDR level
#' @param offset Offset parameter when estimated the FDR. Setting c=1 results in a more conservative FDR result.
#'
#' @return A list containing the vector rej, which are the indices of the columns of X corresponding to the non-null variables,
#' the vector W which are the test statistics for each variable, and the numeric threshold used to create the rejection set.
#'
#' @importFrom knockoff stat.lasso_lambdasmax knockoff.threshold
#' @importFrom SNPknock knockoffDMC
#'
#' @export
knockDMC = function(y, X, q=0.1, offset=0) {
  p = dim(X)[2]

  # First estimate the transition kernel and the initial state distribution.
  # Here we use the MLE
  initial_states = as.numeric(X[, 1])
  initial_states = table(initial_states) + 1
  initial_states = initial_states / sum(initial_states)
  Q = array(dim=c(p-1, 3, 3))
  for (i in 1:(p-1)) {
    subQ = table( c(X[,i]), c(X[,i+1])) + 1
    subQ = diag(1 / rowSums(subQ)) %*% subQ
    Q[i, , ] = subQ
  }

  # knockoffDMC procedure
  Xk = knockoffDMC(X, initial_states, Q)
  W = stat.lasso_lambdasmax(X, Xk, y)
  threshold = knockoff.threshold(W, fdr=q, offset=offset)
  rej = which(W >= threshold)

  return(list(W=W, rej=rej, threshold=threshold))
}





#' Gaussian mirror
#'
#' @description This is a variable selection procedure developed by Xing, Zhao, and Liu (2023)
#' that controls the false discovery rate. This implementation corresponds to Algorithm 1 in their paper.
#' However, instead of using OLS, we use Lasso regression, which boosts power. If you seek to use OLS, just enter lam.lasso=0.
#'
#' @param y Response variable
#' @param X Design matrix
#' @param q The desired FDR level
#' @param lam.lasso The lasso parameter used to fit the regression model. If null, the function will find the optimal
#' value using cross-validation.
#' @param offset Offset parameter when estimated the FDR. Setting c=1 results in a more conservative FDR result.
#'
#' @return A list containing the vector rej, which are the indices of the columns of X corresponding to the non-null variables,
#' the vector W which are the test statistics for each variable, and the numeric threshold used to create the rejection set.
#'
#' @importFrom knockoff knockoff.threshold
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom stats coef rnorm
#'
#' @export
gm_low_d = function(y, X, q=0.1, lam.lasso=NULL, offset=0) {
  # extract dimensionality
  # while OATK normalizes columns of X to unit norm,
  # the GM package normalizes them to \sqrt(n-1) norm,
  # so we will follow the same scaling here.
  n = length(y)
  p = dim(X)[2]
  X = scale(X)

  #  "Low_d GM is only supported for low-dimensional setting (p<n)."
  stopifnot(p < n)

  # compute lasso coefficient, if needed.
  if (is.null(lam.lasso)) {
    # lasso regression on original data
    lambda_max = 1 * max(abs(t(X) %*% y))/(n^(3/2))
    lambda_min = lambda_max/100
    nlambda = 100
    k = (0:(nlambda - 1))/nlambda
    lambda = lambda_max * (lambda_min/lambda_max)^k
    fit0.lasso = cv.glmnet(X, y, family = 'gaussian', lambda = lambda,
      nfolds = 10, alpha = 1)
    lam.lasso = fit0.lasso$lambda.min
  }


  # SVD
  SVD.X = svd(X) #X = U%*%diag(D)%*%t(V)
  V = SVD.X$v #pxp matrix
  D = SVD.X$d #px1 vector
  lam = D^2 #eigenvalues of Gram matrix
  invgram = V%*%((1/lam)*t(V))

  # inialize GM variables by drawing iid gaussians of unit norm
  Z = matrix(rnorm(n * p),n)
  Z = scale(Z) / sqrt(n-1)

  # construct gaussian mirror for each variable and regress
  W = sapply(1:p, function(j){
    invgram.j = invgram[-j,-j] - invgram[j,-j]%*%t(invgram[j,-j])/invgram[j,j] #inv gram w/o x_j, using block matrix inverse result
    XZ.perp = cbind(X[,j], Z[,j]) - X[,-j]%*%(invgram.j%*%(t(X[,-j])%*%cbind(X[,j], Z[,j]))) #Component of [x_j, z_j] orthogonal to X[,-j]
    c.j.num = as.numeric(t(XZ.perp[,1])%*%XZ.perp[,1])
    c.j.den = as.numeric(t(XZ.perp[,2])%*%XZ.perp[,2])
    c.j = sqrt(c.j.num/c.j.den) #constant c_j

    X_GM = cbind(X[, j]+c.j*Z[, j], X[, j]-c.j*Z[, j], X[, -j])

    # lasso
    mdl.lasso = glmnet(X_GM, y, family='gaussian', lambda=lam.lasso, alpha=1)
    beta_hat = coef(mdl.lasso, s='lambda.min')[2:3]
    return(abs(beta_hat[1]+beta_hat[2]) - abs(beta_hat[1]-beta_hat[2]))
  })

  threshold = knockoff.threshold(W, fdr=q, offset=offset)
  rej = which(W >= threshold)

  return(list(W=W, rej=rej, threshold=threshold))
}




#' One-at-a-time Knockoffs (OATK)
#'
#' @description OATK is a variable selection procedure that controls the false discovery rate (FDR).
#' Proposed by Guan, Ren, and Apley (2025), it achieves higher power compared to other FDR-controlling
#' procedures such as Benjamini-Hochberg, knockoff filter, and Gaussian mirror. This implementation conducts
#' a fast version of OATK using ridge regression. If you desire to run OATK with OLS, set lam.ridge=0.
#'
#' @param y Response variable
#' @param X Design matrix
#' @param q The desired FDR level
#' @param lam.ridge The ridge parameter used to fit the regression model. If null, the function will find the optimal
#' value using leave-one-out cross-validation.
#' @param offset Offset parameter when estimated the FDR. Setting c=1 results in a more conservative FDR result.
#' @param U the matrix U from the singular value decomposition of X=UDV^t. If any of U, D, or V is null, the function will compute them
#' using svd().
#' @param D the matrix D from the singular value decomposition of X=UDV^t
#' @param V the matrix V from the singular value decomposition of X=UDV^t
#'
#' @return A list containing the vector rej, which are the indices of the columns of X corresponding to the non-null variables,
#' the vector W which are the test statistics for each variable, the numeric threshold used to create the rejection set,
#' and numeric lam.ridge which was the used ridge parameter.
#'
#'
#' @importFrom knockoff knockoff.threshold
#' @importFrom stats rnorm
#'
#' @export
oatk = function (y, X, q=0.1, lam.ridge=NULL, offset=0, U=NULL, V=NULL, D=NULL)
{
  n = length(y)
  p = dim(X)[2]
  #  "OATK is only supported for low-dimensional setting (p<n)."
  stopifnot(p < n)

  # make sure X is scaled to unit norm
  X = scale(X) / sqrt(n-1)  ## normalize the design matrix
  y = y - mean(y)


  if (is.null(U) || is.null(V) || is.null(D)) {
    # get relevant matrices from SVD
    SVD.X = svd(X) #X = U%*%diag(D)%*%t(V)
    U = SVD.X$u #nxp matrix
    V = SVD.X$v #pxp matrix
    D = SVD.X$d #px1 vector
  }

  lam = D^2 #eigenvalues of Gram matrix

  # "X must be full-rank for OATK!"
  stopifnot(min(lam) > 1e-6)


  # run ridge regression on original X
  # if lam.ridge is NULL, it finds the best lambda
  # using leave-one-out cross-validation
  if (is.null(lam.ridge)) {
    # run ridge
    lam.ridge = max(D^2)*(2/3)^(0:40) #set to 0 for OLS, or a fixed scalar for ridge with specified regularization constant
    mdl.ridge = ridge.reg(Y=y, U=U, D=D, V=V, lam=lam.ridge)
    beta_hat = mdl.ridge$beta.hat  #ridge regression coefficients
    lam.ridge = mdl.ridge$lam.best  #final ridge regression regularization parameter
  } else {
    beta_hat = ridge.reg(Y=y, U=U, D=D, V=V, lam=lam.ridge)$beta.hat
  }

  # construct inverse Gram matrices and OLS coefficient
  Sig.inv.lam =  apply(V, 1, function(x) sum(x^2/(lam+lam.ridge))) #diagonals of the inverse Gram for ridge
  Sig.inv =  apply(V, 1, function(x) sum(x^2/lam)) #diagonals of the inverse Gram for OLS
  beta_OLS = ridge.reg(Y=y, U=U, D=D, V=V, lam=0)$beta.hat #OLS coefficients

  # construct knockoff coefficients
  beta_tilde.fixed = beta_hat - Sig.inv.lam/Sig.inv*beta_OLS #beta_tilde prior to adding the component due to Z
  Z = matrix(rnorm(n * p),n)
  Z = Z - U %*% (t(U) %*% Z) #orthogonalize Z w.r.t. X
  Z = scale(Z) / sqrt(n-1) #scale Z to have unit norm
  r.Z.Y = as.numeric(t(Z)%*%y)
  beta_tilde = as.numeric(beta_tilde.fixed + Sig.inv.lam/sqrt(Sig.inv)*r.Z.Y)

  # generate knockoff statistics and reject
  W = pmax(abs(beta_hat), abs(beta_tilde)) * sign(abs(beta_hat) - abs(beta_tilde))
  threshold = knockoff.threshold(W, fdr=q, offset=offset)
  rej = which(W >= threshold)

  return(list(W=as.numeric(W), rej=rej, threshold=threshold, lam.ridge=lam.ridge))
}



#' One-at-a-time Knockoffs (OATK)
#'
#' @description OATK is a variable selection procedure that controls the false discovery rate (FDR).
#' Proposed by Guan, Ren, and Apley (2025), it achieves higher power compared to other FDR-controlling
#' procedures such as Benjamini-Hochberg, knockoff filter, and Gaussian mirror. This implementation conducts
#' lasso for the knockoff regression. 
#'
#' @param y Response variable
#' @param X Design matrix
#' @param q The desired FDR level
#' @param lam.ridge The ridge parameter used to fit the regression model. If null, the function will find the optimal
#' value using leave-one-out cross-validation.
#' @param offset Offset parameter when estimated the FDR. Setting c=1 results in a more conservative FDR result.
#' @param U the matrix U from the singular value decomposition of X=UDV^t. If any of U, D, or V is null, the function will compute them
#' using svd().
#' @param D the matrix D from the singular value decomposition of X=UDV^t
#' @param V the matrix V from the singular value decomposition of X=UDV^t
#'
#' @return A list containing the vector rej, which are the indices of the columns of X corresponding to the non-null variables,
#' the vector W which are the test statistics for each variable, the numeric threshold used to create the rejection set,
#' and numeric lam.ridge which was the used ridge parameter.
#'
#'
#' @importFrom knockoff knockoff.threshold
#' @importFrom stats rnorm
#'
#' @export
oatk_lasso = function (y, X, q=0.1, lam.ridge=NULL, offset=0, U=NULL, V=NULL, D=NULL)
{
  n = length(y)
  p = dim(X)[2]
  #  "OATK is only supported for low-dimensional setting (p<n)."
  stopifnot(p < n)

  # make sure X is scaled to unit norm
  X = scale(X) / sqrt(n-1)  ## normalize the design matrix
  y = y - mean(y)


  if (is.null(U) || is.null(V) || is.null(D)) {
    # get relevant matrices from SVD
    SVD.X = svd(X) #X = U%*%diag(D)%*%t(V)
    U = SVD.X$u #nxp matrix
    V = SVD.X$v #pxp matrix
    D = SVD.X$d #px1 vector
  }

  lam = D^2 #eigenvalues of Gram matrix

  # "X must be full-rank for OATK!"
  stopifnot(min(lam) > 1e-6)

  # use same lambdas as GM
  lambda_max = 1 * max(abs(t(X) %*% y)) / n
  lambda_min = lambda_max/100
  nlambda = 100
  k = (0:(nlambda - 1))/nlambda
  lambda = lambda_max * (lambda_min/lambda_max)^k
  fit0.lasso = cv.glmnet(X, y, family = 'gaussian', lambda = lambda,
                         nfold = 10, alpha = 1)
  
  lambda_ast = fit0.lasso$lambda.min
  beta_hat = coef(fit0.lasso, s = "lambda.min")[-1]
  
  
  # generate and run knockoffs
  Sig.inv = V %*% diag(1 / lam) %*% t(V)

  sigma = sqrt(1 / diag(Sig.inv))
  B = - Sig.inv %*% diag(sigma^2)
  diag(B) = 0
  proj = X %*% B

  Z <- matrix(rnorm(n * p),n)
  Z <- Z - U %*% (t(U) %*% Z) #orthogonalize Z w.r.t. X
  Z <- scale(Z) / sqrt(n-1) #scale Z to have unit norm
  Z = Z %*% diag(sigma)

  x_tilde = proj + Z
    
  para_list = foreach(j = 1:p, .packages = "glmnet") %do% 
    {
      xnew = cbind(x_tilde[, j], # first column is the knockoff column
                   X[, -j])
      
      # knockoff regression
      fit1.lasso = glmnet(xnew, y, family = 'gaussian', lambda = lambda_ast)
      beta_tilde = coef(fit1.lasso, s = "lambda.min")[2] # accounting for intercept here
      
      return(beta_tilde)
    }
  
  beta_tilde = unlist(para_list)


  # generate knockoff statistics and reject
  W = pmax(abs(beta_hat), abs(beta_tilde)) * sign(abs(beta_hat) - abs(beta_tilde))
  threshold = knockoff.threshold(W, fdr=q, offset=offset)
  rej = which(W >= threshold)

  return(list(W=as.numeric(W), rej=rej, threshold=threshold, lam.ridge=lam.ridge))
}




# OAATK implementation using Ridge
# with a pre-screening procedure using Lasso
# this is the improved version with the rank-2 update of inverse of Sigma_Lam
# with orthogonalization on null knocoffs wrt X_R
oatk_screen_ridge = function (y, x, q = 0.1, offset=0)
{
  n = length(y)
  p = dim(x)[2]

  ############## Pre-screening using Lasso #################
  # lasso regression for screening
  # use same lambdas as GM
  lambda_max = 1 * max(abs(t(x) %*% y)) / n
  lambda_min = lambda_max/100
  nlambda = 100
  k = (0:(nlambda - 1))/nlambda
  lambda = lambda_max * (lambda_min/lambda_max)^k
  fit0.lasso = cv.glmnet(x, y, family = "gaussian", lambda = lambda,
                         nfolds = 10, alpha = 1)


  # obtain the coefficients
  coef_ast = coef(fit0.lasso, s = "lambda.min")[-1]

  # ensure the number of selected variables is less than n
  if (sum(coef_ast != 0) > 0.9 * n) {
    ord_idx = order(abs(coef_ast[which(coef_ast != 0)]))
    coef_ast[which(coef_ast != 0)[1:(sum(coef_ast != 0) -
                                       0.9 * n)]] = 0
  }
  bool.screened = coef_ast != 0
  idx.screened = which(bool.screened)
  num.screened = length(idx.screened)

  if (length(idx.screened)==0 & (p > n)) {
    warning('Pre-screening did not select any variables. Returning null set')
    return(list(W=NULL, rej=NULL, threshold=NULL, lam.ridge=NULL))
  }

  ########### Obtaining the original ridge coefficients ##############
  # run ridge regression and compute inverse ridge gram matrix
  SVD.X <- svd(x) #X = U%*%diag(D)%*%t(V)
  U <- SVD.X$u #nxp matrix
  V <- SVD.X$v #pxp matrix
  D <- SVD.X$d #px1 vector
  lam <- D^2 #eigenvalues of Gram matrix
  lam.ridge <- max(D^2)*(2/3)^(0:40) #set to 0 for OLS, or a fixed scalar for ridge wimatth specified regularization constant


  # high-dimensional --> use 10-fold CV
  if (p > n) {
    fit.ridge = ridge.reg.Kfold(y, x, U, V, D, K=10L, lam.ridge)
    lam.ridge = fit.ridge$lam.best
    beta_hat = fit.ridge$beta.hat
    Sig.inv.lam = qr.solve(t(x)%*%x + lam.ridge * diag(nrow=p) )
    
    # low-dimensional --> use LOOCV
  } else {
    mdl.ridge <- ridge.reg(Y=y, U=U, D=D, V=V, lam = lam.ridge)
    beta_hat <- mdl.ridge$beta.hat  #ridge regression coefficients
    lam.ridge <- mdl.ridge$lam.best  #final ridge regression regularization parameter
    Sig.inv.lam = V %*% diag(1 / (lam+lam.ridge)) %*% t(V)
  }
  beta_hat = as.numeric(beta_hat)

  # Run SVD of the remaining columns
  X_R = x[, bool.screened]
  X_S = x[, !bool.screened]
  SVD.X <- svd(X_R)
  U_R <- SVD.X$u #nxp matrix
  proj_U_R = U_R %*% t(U_R)
  V_R <- SVD.X$v #pxp matrix
  D_R <- SVD.X$d #px1 vector
  lam_R <- D_R^2 #eigenvalues of Gram matrix
  Sig.inv = V_R %*% diag(1 / lam_R) %*% t(V_R)

  ##### Compute knockoff variable for the remaining columns
  # compute the projection of each of the remaining remaining columns to the rest of the remaining column
  sigma_R = sqrt(1 / diag(Sig.inv))
  names(sigma_R) = idx.screened
  B = - Sig.inv %*% diag(sigma_R^2)
  diag(B) = 0
  proj_R = X_R %*% B

  Z_R = matrix(rnorm(n*num.screened), n)
  Z_R <- Z_R - U_R %*% (t(U_R) %*% Z_R) #orthogonalize Z w.r.t. X_R
  Z_R <- scale(Z_R) / sqrt(n-1) #scale Z to have unit norm
  Z_R = Z_R %*% diag(sigma_R)

  x_tilde_R = proj_R + Z_R
  colnames(x_tilde_R) = idx.screened


  ############### Run knockoff regressions ###################
  Xy = t(x) %*% y

  # cl <- makeCluster(ncores)
  # registerDoParallel(cl)
  beta_tilde = foreach(j = 1:p) %do%
    {
      xj = x[, j]
      # generate the knockoff column
      # case 1: screened in variable
      if (j %in% idx.screened) {
        xj_tilde = x_tilde_R[, as.character(j)]
        # case 2: screened out variable --> Need to construct knockoff variable for it
      } else {
        proj = proj_U_R %*% xj
        resid = xj - proj
        sigma_j = l2_norm(resid)

        u = xj - proj_U_R %*% xj
        u = u / l2_norm(u)
        orthog = cbind(U_R, u)

        zj = rnorm(n)
        zj = zj - orthog %*% (t(orthog) %*% zj) 
        zj = zj / sqrt(sum(zj^2))
        xj_tilde = proj + zj * sigma_j
      }

      # Compute inverse of Sigma_inv with the knockoff column
      # note the gram matrix is a rank-2 update, so we use the Woodbury Formula
      ej = numeric(p)
      ej[j] = 1
      v = as.numeric(t(x) %*% (xj_tilde - xj))
      v[j] = 0
      A = cbind(v, ej)
      B = rbind(ej, v)
      K = Sig.inv.lam %*% A
      Sig.inv.lam.tilde = Sig.inv.lam - K %*% solve(diag(nrow=2) + B%*%K, t(K)[2:1, ])

      # calculate knockoff coefficient
      Xy_tilde = Xy
      Xy_tilde[j] = t(xj_tilde) %*% y
      beta_j = as.numeric(Sig.inv.lam.tilde[j, ] %*% Xy_tilde)

      return(beta_j)
    }
  # stopCluster(cl)
  beta_tilde = unlist(beta_tilde)

  W = pmax(abs(beta_hat), abs(beta_tilde)) * sign(abs(beta_hat) - abs(beta_tilde))
  threshold = knockoff.threshold(W, fdr=q, offset=offset)
  rej = which(W >= threshold)


  return(list(W=as.numeric(W), beta_hat=as.numeric(beta_hat), beta_tilde=as.numeric(beta_tilde), rej=rej, threshold=threshold, variable=fifelse(coef_ast != 0, 'remaining', 'screened')))
}





# OAATK implementation using Lasso
# with a pre-screening procedure using Lasso
# In this variant, we do orthogonalize the stochastic component of the null knockoffs
# wrt to X_R.
oatk_screen_lasso = function (y, x, q = 0.1, offset=0)
{
  n = length(y)
  p = dim(x)[2]

  ############## Pre-screening using Lasso #################
  # lasso regression for screening
  # use same lambdas as GM
  lambda_max = 1 * max(abs(t(x) %*% y)) / n
  lambda_min = lambda_max/100
  nlambda = 100
  k = (0:(nlambda - 1))/nlambda
  lambda = lambda_max * (lambda_min/lambda_max)^k
  fit0.lasso = cv.glmnet(x, y, family = "gaussian", lambda = lambda,
                         nfolds = 10, alpha = 1)


  # obtain the coefficients
  beta_hat = coef(fit0.lasso, s = "lambda.min")[-1]
  lam.lasso = fit0.lasso$lambda.min

  # ensure the number of selected variables is less than n
  if (sum(beta_hat != 0) > 0.9 * n) {
    ord_idx = order(abs(beta_hat[which(beta_hat != 0)]))
    beta_hat[which(beta_hat != 0)[1:(sum(beta_hat != 0) -
                                       0.9 * n)]] = 0
  }
  bool.screened = beta_hat != 0
  idx.screened = which(bool.screened)
  num.screened = length(idx.screened)


  if (length(idx.screened)==0 & (p > n)) {
    warning('Pre-screening did not select any variables. Returning null set')
    return(list(W=NULL, rej=NULL, threshold=NULL, lam.ridge=NULL))
  }

  # Run SVD of the remaining columns
  X_R = x[, bool.screened]
  X_S = x[, !bool.screened]
  SVD.X <- svd(X_R)
  U_R <- SVD.X$u #nxp matrix
  proj_U_R = U_R %*% t(U_R)
  V_R <- SVD.X$v #pxp matrix
  D_R <- SVD.X$d #px1 vector
  lam_R <- D_R^2 #eigenvalues of Gram matrix
  Sig.inv = V_R %*% diag(1 / lam_R) %*% t(V_R)

  ##### Compute knockoff variable for the remaining columns
  # compute the projection of each of the remaining remaining columns to the rest of the remaining column
  sigma_R = sqrt(1 / diag(Sig.inv))
  names(sigma_R) = idx.screened
  B = - Sig.inv %*% diag(sigma_R^2)
  diag(B) = 0
  proj_R = X_R %*% B

  Z_R = matrix(rnorm(n*num.screened), n)
  Z_R <- Z_R - U_R %*% (t(U_R) %*% Z_R) #orthogonalize Z w.r.t. X_R
  Z_R <- scale(Z_R) / sqrt(n-1) #scale Z to have unit norm
  Z_R = Z_R %*% diag(sigma_R)

  x_tilde_R = proj_R + Z_R
  colnames(x_tilde_R) = idx.screened


  ############### Run knockoff regressions ###################
  Xy = t(x) %*% y

  # cl <- makeCluster(ncores)
  # registerDoParallel(cl)

  beta_tilde = foreach(j = 1:p, .packages="glmnet") %do%
    {
      xj = x[, j]
      # generate the knockoff column
      # case 1: screened in variable
      if (j %in% idx.screened) {
        xj_tilde = x_tilde_R[, as.character(j)]
      # case 2: screened out variable --> Need to construct knockoff variable for it
      } else {
        proj = proj_U_R %*% xj
        resid = xj - proj
        sigma_j = l2_norm(resid)

        u = xj - proj_U_R %*% xj
        u = u / l2_norm(u)
        orthog = cbind(U_R, u)

        zj = rnorm(n)
        zj = zj - orthog %*% (t(orthog) %*% zj)
        zj = zj / sqrt(sum(zj^2))
        xj_tilde = proj + zj * sigma_j
      }

      x_tilde = cbind(xj_tilde, x[, -j])
      fit1.lasso = glmnet(x_tilde, y, family = "gaussian", lambda = lam.lasso)
      beta = coef(fit1.lasso, s = "lambda.min")[2]


      return(beta)
    }
  # stopCluster(cl)
  beta_tilde = unlist(beta_tilde)

  W = pmax(abs(beta_hat), abs(beta_tilde)) * sign(abs(beta_hat) - abs(beta_tilde))
  threshold = knockoff.threshold(W, fdr=q, offset=offset)
  rej = which(W >= threshold)

  return(list(W=as.numeric(W), beta_hat=as.numeric(beta_hat), beta_tilde=as.numeric(beta_tilde), rej=rej, threshold=threshold, variable=fifelse(beta_hat != 0, 'remaining', 'screened')))

}




#' Derandomized One-at-a-time Knockoffs (OATK)
#'
#' @description This conducts the derandomized version of OATK, which is a variable selection procedure that controls the false discovery rate (FDR)
#' proposed by Guan, Ren, and Apley (2025). They showed derandomization reduces the variances of false discovery proportions and true discovery proportions
#' compared to the basic OATK.
#'
#' @param y Response variable
#' @param X Design matrix
#' @param q The desired FDR level
#' @param lam.ridge The ridge parameter used to fit the regression model. If null, the function will find the optimal
#' value using leave-one-out cross-validation.
#' @param offset Offset parameter when estimated the FDR. Setting c=1 results in a more conservative FDR result.
#' @param M The number of replicates to run for each variable
#' @param eta The frequency threshold to include in the rejection set.
#'
#' @return A list containing the vector rej, which are the indices of the columns of X corresponding to the non-null variables,
#' the matrix rej.mat which contains all the rejection results for each variable across the replicates.
#'
#' @importFrom knockoff knockoff.threshold
#' @importFrom stats rnorm
#'
#' @export
oatk_derandomized = function(y, X, q=0.1, lam.ridge=NULL, M=31, eta=0.5, offset=0)
{
  n = length(y)
  p = dim(X)[2]
  #  "OATK is only supported for low-dimensional setting (p<n)."
  stopifnot(p < n)

  # make sure X is scaled to unit norm
  X = scale(X) / sqrt(n-1)  ## normalize the design matrix
  y = y - mean(y)

  # get relevant matrices from SVD
  SVD.X = svd(X) #X = U%*%diag(D)%*%t(V)
  U = SVD.X$u #nxp matrix
  V = SVD.X$v #pxp matrix
  D = SVD.X$d #px1 vector
  lam = D^2 #eigenvalues of Gram matrix

  # "X must be full-rank for OATK!"
  stopifnot(min(lam) > 1e-6)


  # run ridge regression on original X
  # if lam.ridge is NULL, it finds the best lambda
  # using leave-one-out cross-validation
  if (is.null(lam.ridge)) {
    # run ridge
    lam.ridge = max(D^2)*(2/3)^(0:40) #set to 0 for OLS, or a fixed scalar for ridge with specified regularization constant
    mdl.ridge = ridge.reg(Y=y, U=U, D=D, V=V, lam=lam.ridge)
    beta_hat = mdl.ridge$beta.hat  #ridge regression coefficients
    lam.ridge = mdl.ridge$lam.best  #final ridge regression regularization parameter
  } else {
    beta_hat = ridge.reg(Y=y, U=U, D=D, V=V, lam=lam.ridge)$beta.hat
  }

  # construct inverse Gram matrices and OLS coefficient
  Sig.inv.lam =  apply(V, 1, function(x) sum(x^2/(lam+lam.ridge))) #diagonals of the inverse Gram for ridge
  Sig.inv =  apply(V, 1, function(x) sum(x^2/lam)) #diagonals of the inverse Gram for OLS
  beta_OLS = ridge.reg(Y=y, U=U, D=D, V=V, lam=0)$beta.hat #OLS coefficients

  # construct the determinististic component of beta_tilde,
  # prior to adding the stochastic componenet Z
  beta_tilde.fixed = beta_hat - Sig.inv.lam/Sig.inv*beta_OLS #beta_tilde prior to adding the component due to Z

  # derandomization procedure of the random component of knockoff coefficients
  rej.mat = matrix(0,M,p) #initialize matrix storing selected variables on each replicate
  for (i in 1:M) {
    Z = matrix(rnorm(n * p),n)
    Z = Z - U %*% (t(U) %*% Z) #orthogonalize Z w.r.t. X
    Z = scale(Z) / sqrt(n-1) #scale Z to have unit norm
    r.Z.Y = as.numeric(t(Z)%*%y)
    beta_tilde = as.numeric(beta_tilde.fixed + Sig.inv.lam/sqrt(Sig.inv)*r.Z.Y)

    W = apply(cbind(abs(beta_hat), abs(beta_tilde)), 1, max)*sign(abs(beta_hat) - abs(beta_tilde))
    tau = knockoff.threshold(W, fdr=q, offset=offset)
    rej.mat[i,which(W >= tau)] = 1
  }
  frac.rej = apply(rej.mat,2,mean) #p-length vector of selection fractions across all M replicates
  rej = which(frac.rej > eta) #Selected set, averaging the selection events

  return(list(rej=rej, rej.mat=rej.mat))
}


# Multiple one-at-a-time knockoffs, which
# generates valid p-values and rejects variables based
# on their p-values. For details, see Guan et al. (2024).
# M is the number of generations per variable.

#' Multi-bit One-at-a-time Knockoffs (OATK)
#'
#' @description This conducts the multi-bit version of OATK, which is a variable selection procedure that controls the false discovery rate (FDR) and
#' produces p-values, proposed by Guan, Ren, and Apley (2025).
#'
#' @param y Response variable
#' @param X Design matrix
#' @param q The desired FDR level
#' @param lam.ridge The ridge parameter used to fit the regression model. If null, the function will find the optimal
#' value using leave-one-out cross-validation.
#' @param offset Offset parameter when estimated the FDR. Setting c=1 results in a more conservative FDR result.
#' @param M The number of replicates to run for each variable. Higher M yields finer p-values.
#' @param gamma Rejection hyperparameter.
#'
#' @return A list containing the vector rej, which are the indices of the columns of X corresponding to the non-null variables,
#' the vector p_val which is the p-values for each variable, and the vector stats which is the test statistics for each variable.
#'
#' @importFrom data.table as.data.table between ":="
#' @importFrom stats rnorm
#'
#' @export
oatk_multiple = function (y, X, q=0.1, lam.ridge=NULL, M=3L, gamma=max(1/(M+1), q/(1+q)), offset=0)
{
  n = length(y)
  p = dim(X)[2]
  #  "OATK is only supported for low-dimensional setting (p<n)."
  stopifnot(p < n)
  stopifnot(between(gamma, q/(q+1), 0.5))


  # make sure X is scaled to unit norm
  X = scale(X) / sqrt(n-1)  ## normalize the design matrix
  y = y - mean(y)

  # get relevant matrices from SVD
  SVD.X = svd(X) #X = U%*%diag(D)%*%t(V)
  U = SVD.X$u #nxp matrix
  V = SVD.X$v #pxp matrix
  D = SVD.X$d #px1 vector
  lam = D^2 #eigenvalues of Gram matrix

  # "X must be full-rank for OATK!"
  stopifnot(min(lam) > 1e-6)
  #"We can only generate up to n-p multi-bit knockoffs per variable!"
  stopifnot(M <= n-p)


  # run ridge regression on original X
  # if lam.ridge is NULL, it finds the best lambda
  # using leave-one-out cross-validation
  if (is.null(lam.ridge)) {
    # run ridge
    lam.ridge = max(D^2)*(2/3)^(0:40) #set to 0 for OLS, or a fixed scalar for ridge with specified regularization constant
    mdl.ridge = ridge.reg(Y=y, U=U, D=D, V=V, lam=lam.ridge)
    beta_hat = mdl.ridge$beta.hat  #ridge regression coefficients
    lam.ridge = mdl.ridge$lam.best  #final ridge regression regularization parameter
  } else {
    beta_hat = ridge.reg(Y=y, U=U, D=D, V=V, lam=lam.ridge)$beta.hat
  }

  # construct inverse Gram matrices and OLS coefficient
  Sig.inv.lam = apply(V, 1, function(x) sum(x^2/(lam+lam.ridge))) #diagonals of the inverse Gram for ridge
  Sig.inv = apply(V, 1, function(x) sum(x^2/lam)) #diagonals of the inverse Gram for OLS
  beta_OLS = ridge.reg(Y=y, U=U, D=D, V=V, lam=0)$beta.hat #OLS coefficients

  # construct the determinististic component of beta_tilde,
  # prior to adding the stochastic componenet Z
  beta_tilde.fixed = beta_hat - Sig.inv.lam/Sig.inv*beta_OLS #beta_tilde prior to adding the component due to Z


  # generate M knockoffs per variable and obtain knockoff coefficients
  data = mapply(function(j){
    Z = matrix(rnorm(n * M),n)
    Z = Z - U %*% (t(U) %*% Z) #orthogonalize Z w.r.t. X
    Z = qr.Q(qr(Z)) # make Z an orthonormal matrix
    r.Z.Y = as.numeric(t(Z)%*%y)
    beta_tilde_j = abs( as.numeric(beta_tilde.fixed[j] + Sig.inv.lam[j]/sqrt(Sig.inv[j])*r.Z.Y) )
    beta_j = abs(beta_hat[j])

    pval = (sum(beta_tilde_j >= beta_j) + 1) / (M + 1)
    return(c(pval, max(c(beta_j,beta_tilde_j))))
    }, 1:p)

  # store each variable's test statistic and p-value
  data = cbind(t(data), 1:p)
  colnames(data) = c('pj', 'Mj', 'id')
  data = as.data.table(data)
  p_val = data$pj
  stats = data$Mj

  # sort by Mj and reject
  data = data[order(Mj, decreasing=TRUE)]
  data[, above_thresh:=pj > gamma]
  data[, below_thresh:=!above_thresh]
  data[, test:=(offset + cumsum(above_thresh)) / pmax(1, cumsum(below_thresh))]
  filter_set = which(data$test <= (1-gamma)/gamma*q)
  if (length(filter_set)==0) return(list(rej=integer(0), p_val=p_val, stats=stats))
  k_hat =  max(filter_set)
  rej = data[1:k_hat]
  rej = sort(rej[pj <= gamma]$id)


  return(list(rej=rej, p_val=p_val, stats=stats))
}

# Conditionally calibrated one-at-a-time knockoffs, which
# achieves finite sampel FDR control.
# For details, see Guan et al. (2024).
# M is the number of replicates in the Monte Carlo simulations.

#' Conditionally calibrated one-at-a-time knockoffs (OATK)
#'
#' @description This conducts conditional calibration on the results from the basic OATK procedure,
#' proposed by Guan, Ren, and Apley (2025). T. It achieves finite FDR control.
#' @param y Response variable
#' @param X Design matrix
#' @param q The desired FDR level
#' @param lst.oatk The output list from the basic oatk procedure ran using oatk()
#' @param M The number of replicates to run for each variable during the Monte Carlo simulations
#' @param offset Offset parameter when estimated the FDR. Setting c=1 results in a more conservative FDR result.
#'
#' @return A list containing the vector rej, which are the indices of the columns of X corresponding to the non-null variables,
#' and the vector e_val which is the e-values for the variables.
#'
#' @importFrom stats lm rnorm
#'
#' @export
calibrateOATK = function(y, X, lst.oatk, offset=0, q=0.1, M=300L) {
  n = length(y)
  p = dim(X)[2]
  #  "OATK is only supported for low-dimensional setting (p<n)."
  stopifnot(p < n)

  # make sure X is scaled to unit norm
  X = scale(X) / sqrt(n-1)  ## normalize the design matrix
  y = y - mean(y)

  # extract OATK results
  W = lst.oatk$W
  tau = lst.oatk$threshold
  rej_oatk = lst.oatk$rej
  lam.ridge = lst.oatk$lam.ridge


  # filtering step to reduce the set of variables to run calibration on:
  # for details, see Guan et al. (2024) and Luo et al. (2022)
  # cond_j is the promising set of variables to run the calibration test on

  # collect p-value for the standard t-test
  mdl = lm(Y ~ .-1, data = data.frame("Y" = y, X))
  p_values = summary(mdl)$coefficients[, 4]

  # run the BH procedure at FDR level 4q
  rej_bh = bh(y, X, q*4)$rej

  set1 = intersect(rej_bh, which(p_values <= q / 2))
  filter_thresh = min(sort(abs(W))[p - length(set1)], tau)
  set2 = which(abs(W) >= filter_thresh)
  cond_j = union(union(set1, set2), rej_oatk)
  cond_j = intersect(cond_j, which(W>0))

  e_values = numeric(p)

  # Run the conditional calibration procedure on the set of promising variables
  # if none, then just return the original OATK variables
  # no need to do calibration if filter returns no additional promising variables
  if (length(cond_j)==0) {
    warning('No promising variables found. Returning rejection set of basic OATK')
    return(list(rej=rej))
  } else {
    # compute e-values
    # get relevant matrices from SVD that we use for conditional calibration
    SVD.X = svd(X) #X = U%*%diag(D)%*%t(V)
    U = SVD.X$u #nxp matrix
    V = SVD.X$v #pxp matrix
    D = SVD.X$d #px1 vector
    lam = D^2 #eigenvalues of Gram matrix
    invgram <- V%*%((1/lam)*t(V))

    e_values_filt = sapply(cond_j, function(j) {
      X_nj = X[, -j]
      invgram.j = invgram[-j,-j] - invgram[j,-j]%*%t(invgram[j,-j])/invgram[j,j] #inv gram w/o x_j, using block matrix inverse result
      Pi_nj = X_nj %*% invgram.j %*% t(X_nj) # projection of X_nj
      Pi_nj_perp = diag(n) - Pi_nj # orthogonal projection of X_nj
      py = Pi_nj %*% y
      py = matrix(rep(py, M), n, M)

      V_res = qr.Q(qr(cbind(X_nj, matrix(0,n,n-p+1))))[,p:n]
      u = matrix(rnorm((n-p+1)*M), n-p+1)
      u = scale(u) / sqrt(n-p+1)
      u = l2_norm(Pi_nj_perp %*% y) * V_res %*% u
      y_gen = py + u

      c_hat = abs(tau / W[j])
      phi_vec = sapply(1:M, function(i) {
          # obtain OATK statistics from the simulated y
          lst.resultsGen = oatk(y_gen[, i], X, q, lam.ridge, offset, U, V, D)
          W_repl =  lst.resultsGen$W
          tau_repl = lst.resultsGen$threshold
          W_repl_j = W_repl[j]
          # calculate phi
          phi = (as.numeric(c_hat*W_repl_j >= tau_repl) - as.numeric(W_repl_j <= -tau_repl)) / max(1, sum(W_repl <= -tau_repl))
          return(phi)})
      phi = mean(phi_vec)
      e_j = ifelse(phi > 0, 0, p / max(1, sum(W <= -tau)))

      return(e_j)})
    e_values[cond_j] = e_values_filt
    # Run e-value benjamini-hochberg procedure
    rej = eBH(e_values, q)
    return(list(rej=rej, e_val=e_values))
  }
}

