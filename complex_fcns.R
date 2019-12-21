# Performs an LRT for activation 
# where the first column of X / first entry of beta (beta0) represents a baseline level 
# and the second column of X / second entry of beta (beta1) represents a task-related 
# activation.  The test is Ho: beta1=0 vs. Ha: beta1!=0.
# -yr and yi are the real and imaginary voxel times (for a single voxel)
# -p is the AR order
# -tol is difference between successive log-likelihood values at which the iterative 
#   parameter estimation algorithm stops.
# -max.iter is the maximum number of iterations allowed in the iterative algorithm.

# The function outputs 
# - par: the MLEs under Ha, which is a list with entries:
#   - beta (2 entries)
#   - alpha (the AR coefficients), with p entries
#   - sigma_R^2, the variance of the real time series
#   - sigma_I^2, the variance of the imaginary time series
#   - rho: the correlation between the real and imaginary time series
#   - theta: the constant phase of the complex-valued time series
# - lrt: the value of the LRT statistic
complex.gen.est.lrt <- function(X, yr, yi, p, tol, max.iter)
{
	out1 <- complex.gen.est.ll(X, yr, yi, p, tol, max.iter) #unres
	out0 <- complex.gen.est.ll(cbind(X[,1]), yr, yi, p, tol, max.iter) #res
	lrt <- 2 * (out1$LL - out0$LL)
	list(par=out1$par, lrt=lrt)
}

# Performs the parameter estimation and calculation of the log-likelihood function 
# used above, 
# using separate functions for the temporally independent (p=0) and AR case (p>=1)
complex.gen.est.ll <- function(X, yr, yi, p, tol, max.iter)
{
  if(p==0)
    out <- est.par.ridep.timeindep(X, cbind(yr, yi), tol, max.iter)
  else if(p>=1)
    out <- est.ri.time.dep(X, yr, yi, p, tol, max.iter)
  else
    stop('Inappropriate p')
  out #list with par(list with beta,alpha,sr2,si2,rho,theta)
	  #          LL
}

# Performs the parameter estimation and calculation of the log-likelihood function 
# for the AR case (p>=1)
# This calls the function 'Rwrapper_est_complex_par_ri_temp_dep' from the file 
# complex_Sig=gen.c
est.ri.time.dep <- function(X, yr, yi, p, tol, max.iter)
{
  n <- nrow(X)
  q <- ncol(X)
  out <- .C('Rwrapper_est_complex_par_ri_temp_dep', as.integer(n), as.integer(q), 
            as.integer(p), as.double(yr), as.double(yi), as.double(X), beta=double(q), 
            theta=double(1), sr2=double(1), si2=double(1), rho=double(1), 
            alpha=double(p), as.double(tol), as.integer(max.iter), 
            LL=double(1))
  if(out$beta[1] < 0){
    out$beta <- -out$beta
    out$theta <- ifelse(out$theta > 0, out$theta - pi, out$theta + pi)
  }
  par <- list(beta=out$beta, alpha=out$alpha, sr2=out$sr2, si2=out$si2, 
	rho=out$rho, theta=out$theta)
  list(par=par, LL=out$LL)
}

# Performs the parameter estimation and calculation of the log-likelihood function 
# for the temporally independent case (p=0)
est.par.ridep.timeindep <- function(X, y, tol, max.iter) #time indep
{
  conv <- 0; iter <- 0
  n <- nrow(X)
  yr <- y[,1]; yi <- y[,2]
  XpX <- t(X) %*% X
  br <- solve(XpX) %*% t(X) %*% yr
  bi <- solve(XpX) %*% t(X) %*% yi
  Brr <- as.numeric(t(br) %*% XpX %*% br)
  Bii <- as.numeric(t(bi) %*% XpX %*% bi)
  Bri <- as.numeric(t(br) %*% XpX %*% bi)
  theta <- 0.5 * atan2(2*Bri, Brr - Bii)
  beta <- br * cos(theta) + bi * sin(theta)
  sr2 <- mean((yr - X %*% beta * cos(theta))^2)
  si2 <- mean((yi - X %*% beta * sin(theta))^2)
  rho <- mean((yr - X %*% beta * cos(theta)) * (yi - X %*% beta * sin(theta))) / sqrt(sr2*si2)
  LL.new <- comp.LL(n, yr, yi, X, beta, theta, sr2, si2, rho)
  while(conv == 0 & iter <= max.iter){
    LL.old <- LL.new
    iter <- iter + 1
    D <- Brr/sr2^2 + rho^2/(sr2*si2)*Bii - 2*rho/(sr2^1.5*si2^.5)*Bri
    E <- Bii/si2^2 + rho^2/(sr2*si2)*Brr - 2*rho/(sr2^.5*si2^1.5)*Bri
    F <- Bri*(1+rho^2)/(sr2*si2) - rho/sqrt(sr2*si2) * (Brr/sr2 + Bii/si2)
    A <- D/si2 - E/sr2
    B <- -F * (1/sr2 + 1/si2) - rho/sqrt(sr2*si2)*(D+E)
    C <- rho/sqrt(sr2*si2) * (D-E) + F*(1/sr2 - 1/si2)
    psi <- atan2(B, A)
    theta <- 0.5 * (asin(C/sqrt(A^2+B^2)) - psi)
    beta <- (br * (cos(theta)/sr2 - rho/sqrt(sr2*si2)*sin(theta)) + 
               bi * (sin(theta)/si2 - rho/sqrt(sr2*si2)*cos(theta))) / 
      (cos(theta)^2/sr2 - 2*rho/sqrt(sr2*si2)*sin(theta)*cos(theta) + sin(theta)^2/si2)
    sr2 <- mean((yr - X %*% beta * cos(theta))^2)
    si2 <- mean((yi - X %*% beta * sin(theta))^2)
    rho <- mean((yr - X %*% beta * cos(theta)) * (yi - X %*% beta * sin(theta))) / sqrt(sr2*si2)
    LL.new <- comp.LL(n, yr, yi, X, beta, theta, sr2, si2, rho)
    if(LL.new - LL.old < tol) conv <- 1
  }
  if(iter > max.iter) print('Warning: Over max.iter')
  par <- list(beta=beta, theta=theta, sr2=sr2, si2=si2, rho=rho)
  list(par=par, LL=LL.new)
}

# Helper function for the above function, which calculates LL function for the p=0 case
comp.LL <- function(n, yr, yi, X, beta, theta, sr2, si2, rho)
{
  -n/2 * log(sr2*si2*(1-rho^2)) - n
}



# Performs the parameter estimation and calculation of the log-likelihood function 
# under the model that assumes the real and imaginary time series are independent with 
# the same variance (i.e. sigma_R^2 = sigma_I^2 = sig2 and rho=0).  
# Either temporal independence (p=0) or AR dependence (p>=1) can be used with this function

# This function calls the function "Rwrapper_complex_activation" 
# from the file 'complex_Sig=sig2I.c'
complex.ri.indep <- function(X, yR, yI, p, max.iter, LL.eps)
{
  n <- nrow(X)
  q <- ncol(X)
  C <- rbind(c(0,1))
  m <- nrow(C)
  len <- q + 2+ p
  out <- .C("Rwrapper_complex_activation", as.integer(n), as.integer(q), 
            as.integer(p), as.integer(m), as.double(as.vector(C)), 
            as.double(as.vector(X)), as.double(yR), as.double(yI), 
            as.integer(max.iter), as.double(LL.eps), par=double(2*len), 
            lrt=double(1))
  par <- array(out$par, dim=c(2, len))
  beta <- par[1,1:2]
  sig2 <- par[1,4]
  theta <- par[1,3]
  if(p==0) alpha <- NA
  else alpha <- par[1,5:(4+p)]
  list(lrt=out$lrt, 
       par=list(beta=beta, alpha=alpha, sig2=sig2, theta=theta))
}

# Simulates a complex-valued time series with supplied X matrix and parameter values
# Returns a matrix with n rows (the same number of rows as X) 
# and 2 columns: column 1 is the real time series and column 2 is the imaginary time series
sim.complex.ts <- function(X, beta, theta, sr, si, rho, alpha)
{
  if(ncol(X) != length(beta)) stop('incompatible X, beta')
  n <- nrow(X)
  y <- array(dim=c(n, 2))
  mu <- X %*% beta
  y[,2] <- mu * sin(theta) + arima.sim(list(ar=alpha), n=n, sd=si)
  yr.mean <- mu * cos(theta) + rho * sr/si * (y[,2] - mu*sin(theta))
  y[,1] <- yr.mean + arima.sim(list(ar=alpha), n=n, sd=sr * sqrt(1-rho^2))
  y
}