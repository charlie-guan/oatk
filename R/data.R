library(stringr)



# generate y given deisgn matrix X, indices and sign of nonnull variables,
# and signal amplitude

#' getY
#'
#' @description Generates the response vector Y based on the linear model.
#'
#' @param X Design matrix
#' @param beta_true The true coefficients of the model
#' @param amp Signal amplitude
#'
#' @return Numeric vector
#' @export
#' @importFrom stats rnorm
getY = function(X, beta_true, amp=1) {
	y = X %*% (beta_true*amp) + rnorm(nrow(X))
	y = y - mean(y)
	return(y)
}

#' getGaussianX
#'
#' @description generates a random design matrix X from a multivariate normal with
#' zero mean and a covariance structure specified by the input arguments. See
#' Guan, Ren and Apley (2025) for further details.
#'
#' @param n Number of rows of X
#' @param p Number of columns of X
#' @param rho Scaling parameter of the covariance matrix
#' @param mode Covariance structure type. Takes either 'power_decay', 'const_pos', or 'const_neg'.
#' 'power_decay' is the power decay covariance with decay rate rho.
#' 'const_pos' is constant positive covariance with constant rho.
#' 'const_neg' is constant negative covariance with constant rho in the precision matrix.
#'
#' @return Design matrix X with n rows and p columns
#' @export
#' @importFrom MASS mvrnorm
getGaussianX = function(n, p, rho=0.4, mode="power_decay") {
	stopifnot(mode %in% c('power_decay', 'const_pos', 'const_neg'))

	if (mode=='power_decay') {
		Sigma = rho ^ outer(1:p, 1:p, function(i, j) abs(i - j))
	} else if (mode=='const_pos') {
		Sigma = matrix(rho, nrow=p, ncol=p)
  		diag(Sigma) = 1
	} else { # mode=='const_neg'
		Q = matrix(rho, nrow=p, ncol=p)
		diag(Q) = 1
		Sigma = solve(Q)
		Sigma = Sigma / diag(Sigma)[1] # make sure diagonals are one
	}

	X = mvrnorm(n, mu=rep(0, p), Sigma)
	X = scale(X) / sqrt(n-1)  ## normalize the design matrix

  	return(X)
}



#' getDMCX
#'
#' @description generates a random design matrix X where each row is an independent sample
#' from a discrete-time Markov Chain model. For details, see Guan, Ren, and Apley (2025)
#'
#' @param n Number of rows of X
#' @param p Number of columns of X
#'
#' @return Design matrix X with n rows and p columns
#' @export
#' @importFrom stats runif
getDMCX = function(n, p) {
	states = 0:2
	gamma = runif(p-1, 0, 0.5)
	X = mapply(function(i) {
		x = numeric(p)
		same_state = runif(p-1) <= 1 / 3 + 2 * gamma / 3
		x[1] = sample(states, size=1)
		for (i in 2:p) {
	    if (same_state[i-1]) {
	        x[i] = x[i-1]
	      } else {
	        x[i] = sample(setdiff(states, x[i-1]), size=1)
	      }
	  	}
	  	return(as.integer(x))
		}, 1:n)
	return(t(X))
}


#' getHIVData
#'
#' @description Creates a design matrix X from Stanford University's
#' HIV Drug Resistance Database. Adopted from https://hivdb.stanford.edu/download/GenoPhenoDatasets/DRMcv.R.
#'
#' @param dataset string indicating which dataset to pull data from. Options are 'PI', 'NRTI', and 'NNRTI'.
#' @param min.muts numeric indicating how many times a mutation must occur in the dataset to be included in X
#'
#' @return Design matrix X
#' @export
#' @importFrom utils read.table
#' @importFrom stringr str_split str_detect str_replace
getHIVData = function(dataset='PI', min.muts=10) {
  ## automatically read in the data using url
  if(dataset=="PI"){
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt",
                      header=TRUE, sep="\t",
                      stringsAsFactors=FALSE)
    dat[dat=="."] <- NA
    posu <- dat[,10:108]
  }
  if(dataset=="NRTI"){
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/NRTI_DataSet.txt",
                      header=TRUE, sep="\t", comment.char="@",
                      stringsAsFactors=FALSE)
    dat[dat=="."] <- NA
    posu <- dat[,8:247]
  }
  if(dataset=="NNRTI"){
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt",
                      header=TRUE, sep="\t", comment.char="@",
                      stringsAsFactors=FALSE)
    dat[dat=="."] <- NA
    posu <- dat[,6:245]
  }

  # manual list of all unique mutations present in the data
  muts.in = unique(unlist(str_split(dat[, "CompMutList"], ', ')))
  muts.in = str_replace(muts.in, '.', '')
  muts.in = muts.in[nchar(muts.in) == 3]
  muts.in = muts.in[!str_detect(muts.in, '\\*')]
  muts.in = muts.in[str_detect(muts.in, "[0-9]{2}")]

  # validate mutations
  muts.in.conv <- convert.muts(muts.in)
  check.muts(muts.in.conv)

  # get the amino acids and positions for the mutations to be included in the model
  mut <- ifelse(nchar(muts.in)==3,toupper(substr(muts.in,3,3)),
                toupper(substr(muts.in,4,4)))
  ps <- suppressWarnings(ifelse(nchar(muts.in)==3,as.numeric(substr(muts.in,1,2)),
                                as.numeric(substr(muts.in,1,3))))


  # construct design matrix for OLS
  X <- buildX(posu, mut, ps)


  # remove all rows with missing values
  rem.rows <- unique(which(is.na(X),arr.ind=TRUE)[,1])
  df.log.cc <- X[-rem.rows,]  # complete case

  # remove mutations that are rare
  rare.muts <- which(colSums(df.log.cc)<min.muts)
  if(length(rare.muts)>0){
    message(paste0(muts.in[rare.muts],
                   " excluded from the model because it appears in fewer than ",
                   min.muts," sequences.\n"))
    df.log.cc <- df.log.cc[,-(rare.muts)]
  }

  n <- nrow(df.log.cc) # number of samples
  df.log.cc = scale(df.log.cc) / sqrt(n-1)
  return(df.log.cc)
}



#' getHIVDataWithY
#'
#' @description Creates a matrix consisting of the response vector y corresponding
#' to the fold resistance data and a design matrix X from Stanford University's
#' HIV Drug Resistance Database. Adopted from https://hivdb.stanford.edu/download/GenoPhenoDatasets/DRMcv.R.
#'
#' @param dataset string indicating which dataset to pull data from. Options are 'PI', 'NRTI', and 'NNRTI'.
#' @param min.muts numeric indicating how many times a mutation must occur in the dataset to be included in X
#' @param drug string indicating which fold resistance data to collect for y.
#' @return list with entry y (the response variable) and X (the design matrix)
#' @export
#' @importFrom utils read.table
#' @importFrom stringr str_split str_detect str_replace
getHIVDataWithY = function(dataset='PI', min.muts=10, drug="FPV") {

  ## automatically read in the data using url
  if(dataset=="PI"){
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt",
                      header=TRUE, sep="\t",
                      stringsAsFactors=FALSE)
    dat[dat=="."] <- NA
    posu <- dat[,10:108]
  }
  if(dataset=="NRTI"){
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/NRTI_DataSet.txt",
                      header=TRUE, sep="\t", comment.char="@",
                      stringsAsFactors=FALSE)
    dat[dat=="."] <- NA
    posu <- dat[,8:247]
  }
  if(dataset=="NNRTI"){
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt",
                      header=TRUE, sep="\t", comment.char="@",
                      stringsAsFactors=FALSE)
    dat[dat=="."] <- NA
    posu <- dat[,6:245]
  }

  # manual list of all unique mutations present in the data
  muts.in = unique(unlist(str_split(dat[, "CompMutList"], ', ')))
  muts.in = str_replace(muts.in, '.', '')
  muts.in = muts.in[nchar(muts.in) == 3]
  muts.in = muts.in[!str_detect(muts.in, '\\*')]
  muts.in = muts.in[str_detect(muts.in, "[0-9]{2}")]

  # validate mutations
  muts.in.conv <- convert.muts(muts.in)
  check.muts(muts.in.conv)

  # get the amino acids and positions for the mutations to be included in the model
  mut <- ifelse(nchar(muts.in)==3,toupper(substr(muts.in,3,3)),
                toupper(substr(muts.in,4,4)))
  ps <- suppressWarnings(ifelse(nchar(muts.in)==3,as.numeric(substr(muts.in,1,2)),
                                 as.numeric(substr(muts.in,1,3))))


  # construct design matrix
  X <- buildX(posu, mut, ps)

  # construct dependent variable
  drugcol <- which(names(dat)==drug)
  Y <- as.numeric(dat[,drugcol])  # absolute measure
  Ylog10 <- log10(Y)
  df.log <- data.frame(Y=Ylog10, X=X)


  # remove all rows with missing values
  rem.rows <- unique(which(is.na(df.log),arr.ind=TRUE)[,1])
  df.log.cc <- df.log[-rem.rows,]  # complete case

  # remove mutations that are rare
  rare.muts <- which(colSums(df.log.cc[, -1])<min.muts) #-1 removes the y column from consideration
  if(length(rare.muts)>0){
    message(paste0(muts.in[rare.muts],
                   " excluded from the model because it appears in fewer than ",
                   min.muts," sequences.\n"))
    df.log.cc <- df.log.cc[,-(rare.muts+1)]
  }

  df.log.cc = df.log.cc[!is.infinite(df.log.cc$Y), ]

  # center and scale data
  y = df.log.cc$Y
  y = y - mean(y)
  n = length(y)
  X = df.log.cc[, -1]
  X = scale(X) / sqrt(n-1)
  return(list(y=y, X=X))
}



#' sampleNonnull
#'
#' @description Creates a random coefficient vector with k nonzero entries, whose signs
#' are sampled uniformly at random.
#'
#' @param p Numeric indicating length of vector
#' @param k Numeric indicating number of nonzero entries
#' @return A vector
#' @export
#' @importFrom stats rbinom
sampleNonnull = function(p, k) {
	stopifnot(k < p)
	beta = rep(0, p)
	beta[sample(1:p, k)] = rbinom(k, 1, 0.5)*2-1
	return(beta)
}

#' calculateFDP
#'
#' @description Calculates the false discovery proportion.
#'
#' @param beta_true The true regression coefficients.
#' @param rej The rejection set
#' @return A numeric indicating the FDP.
#' @export
calculateFDP = function(beta_true, rej) {
	if (is.null(rej)) return(0) else return( sum( beta_true[rej]==0 ) / max(1, length(rej)) )
}

#' calculateTDP
#'
#' @description Calculates the true discovery proportion.
#'
#' @param beta_true The true regression coefficients.
#' @param rej The rejection set
#' @return A numeric indicating the TDP
#' @export
calculateTDP = function(beta_true, rej) {
	if (is.null(rej)) return(0)
	k = sum(beta_true!=0)
	return( sum( beta_true[rej]!=0 ) / k )
}


### Helper function to load the HIV Database
#' convert.muts
#'
#' @description Helper function to check that the mutations have been entered correctly
#' if the letters are entered as lower case, they are converted to upper case. Convert insertions or deletions to # or ~.
#' Taken from https://hivdb.stanford.edu/download/GenoPhenoDatasets/DRMcv.R
#' @param muts.in Input mutation labels
#'
#' @return Corrected input mutation labels
convert.muts <- function(muts.in){
  muts.in5 <- which(nchar(muts.in) == 5)
  muts.in6 <- which(nchar(muts.in) == 6)
  for(mi in muts.in5){
    postmp <- substr(muts.in[mi],1,2)
    if(substr(muts.in[mi],3,5) == "ins") muts.in[mi] <- paste0(postmp,"#")
    if(substr(muts.in[mi],3,5) == "del") muts.in[mi] <- paste0(postmp,"~")
  }
  for(mi in muts.in6){
    postmp <- substr(muts.in[mi],1,3)
    if(substr(muts.in[mi],4,6) == "ins") muts.in[mi] <- paste0(postmp,"#")
    if(substr(muts.in[mi],4,6) == "del") muts.in[mi] <- paste0(postmp,"~")
  }
  return(muts.in)
}

#' check.muts
#'
#' @description Helper function to check that the mutations have been entered correctly.
#' Taken from https://hivdb.stanford.edu/download/GenoPhenoDatasets/DRMcv.R
#' @param muts.in Input mutation labels
#'
#' @return Nothing if all mutations are correctly labeled. Otherwise, it throws an error.
check.muts <- function(muts.in){
  # all entries should be nchar=3
  if(any(nchar(muts.in) < 3) | any(nchar(muts.in) > 4))
    stop('All entries in argument "muts.in" should be between 3 and 4 characters long.')

  muts.in3 <- muts.in[nchar(muts.in) == 3]
  muts.in4 <- muts.in[nchar(muts.in) == 4]

  # all should have numbers first two or three characters and a letter for the last
  if(!all(toupper(substr(muts.in3,3,3))%in%c(LETTERS,"#","~")))
    stop('All entries in argument "muts.in" must have a letter, #, or ~ in the last
         character.')
  if(any(is.na(as.numeric(substr(muts.in3,1,2)))))
    stop('All entries in argument "muts.in" must begin in two or three digits.')

  if(!all(toupper(substr(muts.in4,4,4))%in%c(LETTERS,"#","~")))
    stop('All entries in argument "muts.in" must have a letter, #, or ~ in the last
         character.')
  if(any(is.na(as.numeric(substr(muts.in4,1,3)))))
    stop('All entries in argument "muts.in" must begin in two or three digits.')
}


#' buildX
#'
#' @description Helper function to create the design matrix X with the input mutations/positions.
#' Taken from https://hivdb.stanford.edu/download/GenoPhenoDatasets/DRMcv.R
#' @param dat Raw data from the HIV Stanfrod Database
#' @param mut list of mutation to keep
#' @param ps position of the mutations
#'
#' @return Design matrix X
buildX <- function(dat, mut, ps){
  X <- matrix(NA, nrow=nrow(dat), ncol=length(mut))

  # loop through all positions
  for(p in unique(ps)){
    p1 <- substr(dat[,p],1,1)  # first mutation at this position
    p2 <- substr(dat[,p],2,2)
    for(ind in which(ps==p)){
      X[,ind] <-  as.numeric(p1==as.character(mut[ind]) |
                               p2==as.character(mut[ind]))
    }
  }
  colnames(X) <- paste0(ps,mut)
  return(X)
}

