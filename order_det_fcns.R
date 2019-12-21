# These functions perform order detection for a complex-valued time series 
# via the CV models and the sequential testing approach described in Section 3.5.
# Functions from "complex_fcns.R" are called to fit the models for particular 
# AR orders.

# This function detects the AR order under the model with the general structure 
# for the covariance of the real/imaginary errors.

# The inputs are:
# - X: design matrix for the expected BOLD magnitude response
# - yr/yi: real/imaginary time series
# - max.iter/tol: settings for convergence of the iterative algorithm for parameter 
#   estimation (also described in complex_fcns.R)
# - pmax: maximum AR order allowed 
# - signif: significance level applied to each of the sequential hypothesis tests
  
# Output: The detected AR order  
order.det.lrt.complex.gen <- function(X, yr, yi, max.iter, tol, pmax, signif)
{
  phat <- 0; det <- 0
  thresh <- qchisq(1-signif, df=1)
  LL0 <- est.par.ridep.timeindep(X, cbind(yr, yi), tol, max.iter)$LL
  k <- 1
  while(det==0 && k<=pmax){
    LL1 <- est.ri.time.dep(X, yr, yi, k, tol, max.iter)$LL
    lrt <- 2 * (LL1 - LL0)
    if(lrt < thresh){
      det <- 1
      phat <- k - 1
    }
    else{
      k <- k+1
      LL0 <- LL1
    }
  }
  if(k>pmax) phat <- pmax
  phat
}



# Similar function that detects the AR order under the model that assumes 
# the real and imaginary time series are independent with 
# the same variance (i.e. sigma_R^2 = sigma_I^2 = sig2 and rho=0).

order.det.lrt.complex.sig2I <- function(X, yR, yI, max.iter, LL.eps, 
	pmax, signif)
{
	n <- nrow(X)
	phat <- 0
	det <- 0
	thresh <- qchisq(1-signif, df=1)
	LL0 <- compute.LL.complex(X, yR, yI, max.iter, LL.eps, phat)
	k <- 1
	while(det==0 && k<=pmax){
		LL1 <- compute.LL.complex(X, yR, yI, max.iter, LL.eps, k)
		stat <- 2 * (LL1 - LL0)
		if(stat < thresh){
			det <- 1
			phat <- k-1
		}
		else{
			k <- k + 1
			LL0 <- LL1
		}
	}
	if(k>pmax) phat <- pmax
	phat
}

# Helper function for the above function that calculates the log-likelihood
# Calls C function "Rwrapper_complex_unres_only" from complex_Sig=sig2I.c
compute.LL.complex <- function(X, yR, yI, max.iter, LL.eps, p)
{
	n <- nrow(X)
	q <- ncol(X)
	len <- q + 2+ p
	out <- .C("Rwrapper_complex_unres_only", as.integer(n), as.integer(q), 
		as.integer(p), as.double(as.vector(X)), as.double(yR), 
		as.double(yI), as.integer(max.iter), as.double(LL.eps), 
		par=double(len), LL=double(1))
	out$LL
}