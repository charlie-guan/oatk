# OATK: One-at-a-time Knockoffs for Controlled Variable Selection with Higher Power


## Installation
```R
install_github("charlie-guan/oatk")
```

## Examples
The following examples illustrate the basic usage of the OATK package. OATK is designed for any design matrix that is low-dimensional, i.e., p (# of columns) < n (# of rows). These examples are inspired by the numerical results of "One-at-a-time knockoffs: controlled false discovery rate with higher power" by Guan, Ren, and Apley (2025). The preprint is available at https://arxiv.org/abs/2502.18750.

## Example 1
We generate the design matrix with i.i.d. rows, and each row is multivariate normal with zero mean and power-decaying covariance. Below, we conduct the basic OATK procedure, derandomized OATK (which reduces variance and improves FDR control), and multi-bit OATK (which generates p-values).


```R
n = 1000  # Number of observations
p = 300  # Number of predictors 
k = 30 # Number of true nonnull predictors
amp = 5 # Signal Amplitude
beta = sampleNonnull(p, k)

# Generate the variables and response
X = getGaussianX(n, p, mode='power_decay')
y = getY(X, beta, amp)

# Conduct (basic) OATK with a desired FDR of 10%.
lst.oatk = oatk(y, X, q=0.1)
rej = lst.oatk$rej

# Conduct derandomized OATK
rej_derand = oatk_derandomized(y, X, q=0.1)$rej

# Conduct multi-bit OATK with 20 replicates. 
rej_mult = oatk_multiple(y, X, q=0.1, M=20)$rej

```



## Example 2
In this example, we generate the design matrix from a Markov chain and conduct OATK. We also conduct conditional calibration on the OATK output, which achieves finite FDR control.
```R
n = 1000  # Number of observations
p = 100  # Number of predictors 
k = 30 # Number of true nonnull predictors
amp = 5 # Signal Amplitude
beta = sampleNonnull(p, k)

# Generate the variables and response
X = getDMCX(n, p)
y = getY(X, beta, amp)

# Conduct (basic) OATK with a desired FDR of 10%.
lst.oatk = oatk(y, X, q=0.1)
rej = lst.oatk$rej

# Conduct conditional calibration on OATK
# Warning: Speed could be vastly improved by parallelizing the procedure.
rej_calibrated =  calibrateOATK(y, X, lst.oatk, q=0.1)$rej
```


## Example 3
In this example, we consider real X and y from Stanford University's HIV Database to conduct a genome-wide association study (GWAS). We benchmark OATK with competing algorithms to identify mutations that impact fold resistance of the drug fosamprenavir. The results correspond to Table 1 of Guan, Ren, and Apley (2025). 

```R
# Load HIV data and set up params
set.seed(1000)
lst = getHIVDataWithY()
y = lst$y
X = lst$X
mutations = colnames(X)
mutations = stringr::str_replace(mutations, 'X.', '')


# variable selection
rej_bh = sort(mutations[bh(y, X)$rej]) # Benjamini-Hochberg
rej_bc = sort(mutations[bc(y, X)$rej]) # Knockoff filter
rej_gm = sort(mutations[gm_low_d(y, X)$rej]) # Gaussian mirror
rej_oatk = sort(mutations[oatk(y, X)$rej]) # One-at-a-tome knockoffs

# print results
length(rej_oatk)
length(rej_gm)
length(rej_bc)
length(rej_bh)

# uniquely identified mutations by each method
paste0(sort(setdiff(rej_oatk, unique(c(rej_gm, rej_bc, rej_bh)))), collapse=', ')
paste0(sort(setdiff(rej_gm, unique(c(rej_oatk, rej_bc, rej_bh)))), collapse=', ')
paste0(sort(setdiff(rej_bc, unique(c(rej_oatk, rej_gm, rej_bh)))), collapse=', ')
paste0(sort(setdiff(rej_bh, unique(c(rej_oatk, rej_gm, rej_bc)))), collapse=', ')
```
