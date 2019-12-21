

# This function performs the removal of trends in the magnitude and phase 
# for a complex-valued time series via the CV running line as described in 
# Section 3.1 of the paper.
# The inputted arguments are:
# - r: the magnitude time series
# - phi: the phase time series
# - tol: the difference between log-likelihood functions for successive iterations 
#   that constitutes convergence (and stops the algorithm)
# - k: the number of nearest neighbors used when determining the neighborhood in 
#   which the running line is fit
# - max.iter: the maximum number of iterations performed
# - ma.line: takes the value 'ma' or 'line'
#   for 'ma', a moving average is used as a smoother (i.e. assuming constant 
#     magnitude and phase)
#   for 'line', the complex-valued running line is used
# - central.phase: the value theta_0 used to determine the detrended phase 
#   time series \tilde{phi}
# - interp: frequency in which the complex-valued running line is estimated.
#   For instance, interp=1 for every time point, interp=2 for at every other 
#   time point, interp=10 for at every 10 time points, etc. (where linear 
#   interpolation determines the smoothed values in between)

# Outputted arguments: list with two values
# - detrend: matrix with n rows and 2 columns containing the detrended values 
#   of the real (column 1) and imaginary (col 2) time series
# - smoothed: matrix with n rows and 2 columns containing the smoothed values 
#   (i.e. fits of CV running line) of the magnitude (col 1) and phase (col 2)
#   time series.
  
cv.remove.trend.ts <- function(r, phi, tol, k, max.iter, ma.line, 
	central.phase, interp)
{
  out <- smooth.ts(r, phi, tol, k, max.iter, ma.line, interp)
  theta <- out$theta
  theta <- ifelse(theta > -pi & theta <pi, theta, 
                  ifelse(theta <= -pi, theta + 2*pi, 
                         ifelse(theta >= pi, theta - 2*pi, NA)))
  r.fil <- r - out$rho + mean(out$rho)
  phi.fil <- phi - out$theta + central.phase
  phi.fil <- ifelse(phi.fil > -pi & phi.fil <pi, phi.fil, 
                    ifelse(phi.fil <= -pi, phi.fil + 2*pi, 
                           ifelse(phi.fil >= pi, phi.fil - 2*pi, NA)))
  list(detrend=cbind(r.fil*cos(phi.fil), r.fil*sin(phi.fil)), 
       smoothed=cbind(out$rho, theta))
}

# Helper function for performing the CV running line
# Calls the C function 'Rwrapper_complex_running_line' in the 
# file complex_smoother.c
smooth.ts <- function(r, phi, tol, k, max.iter, ma.line, interp)
{
  N <- length(r)
  if(ma.line=='ma') type <- 0
  if(ma.line=='line') type <- 1
  out <- .C('Rwrapper_complex_running_line', as.integer(N), as.integer(k), 
            as.double(r), as.double(phi), as.double(tol), rho=double(N), 
            theta=double(N), n.iter=integer(N), as.integer(max.iter), 
            as.integer(type), as.integer(interp))
  rho <- out$rho; theta <- out$theta
  #interpolation
  if(interp >=2){
	mu.r <- rho * cos(theta)
	mu.i <- rho * sin(theta)
	t <- 1:N
	rem <- (t-1) %% interp
	floor <- interp * ((t-1) %/% interp) + 1
	non.interp <- (N-interp+2):N
	mur.interp <- c(ifelse(rem==0, mu.r, (1-rem/interp)*mu.r[floor] + 
                       rem/interp*mu.r[floor+interp])[-non.interp], 
					   mu.r[non.interp])
	mui.interp <- c(ifelse(rem==0, mu.i, (1-rem/interp)*mu.i[floor] + 
                       rem/interp*mu.i[floor+interp])[-non.interp], 
					   mu.i[non.interp])
	rho <- sqrt(mur.interp^2 + mui.interp^2)
	theta <- atan2(mui.interp, mur.interp)
  }
  list(rho=rho, theta=theta)
}



# Performs similar detrending via a running line for a magnitude-only 
# time series.
# Calls the C function 'Rwrapper_mag_only_run_line' in the 
# file complex_smoother.c

# Output: list with 2 objects:
# - rho: fitted running line
# - detrend: detrended magnitude time series
mag.running.line <- function(r, k, interp)
{
  N <- length(r)
  rho <- .C('Rwrapper_mag_only_run_line', as.integer(N), as.integer(k), 
            as.double(r), rho=double(N), as.integer(interp))$rho
  #interpolation
  if(interp >=2){
	t <- 1:N
	rem <- (t-1) %% interp
	floor <- interp * ((t-1) %/% interp) + 1
	non.interp <- (N-interp+2):N
	rho.interp <- c(ifelse(rem==0, rho, (1-rem/interp)*rho[floor] + 
                       rem/interp*rho[floor+interp])[-non.interp], 
					   rho[non.interp])
	rho <- rho.interp
  }
  list(rho=rho, detrend=r - rho + mean(rho))
}




# Used for spatial smoothing


# Calculates the 3-dimensional array of weights for a 3-dimensional 
# discrete isotropic Gaussian filter from the fwhm in voxels
# Output: normalized 3-dimensional array
kern.3d <- function(fwhm)
{
  #requires package spatialfil
  sig <- fwhm / (2*sqrt(2*log(2)))
  kerndim <- nrow(convKernel(sigma = sig, k = "gaussian")$matrix)
  kw <- .5 * (kerndim - 1)
  kw2 <- -kw:kw
  kern <- array(dim=c(kerndim, kerndim, kerndim))
  for(i in 1:kerndim) for(j in 1:kerndim) for(k in 1:kerndim)
    kern[i,j,k] <- exp(-1/(2*sig^2)*(kw2[i]^2 + kw2[j]^2 + kw2[k]^2))
  kern.norm <- kern / sum(kern)
  kern.norm
}

# Performs the spatial smoothing of a 3-dimensional image (at a single 
# time point).
# Calls the C function 'Rwrapper_spatial_smooth3d' from 
# the file spatial_smooth.c
# Inputs: 
# - image: 3-dimensional array of data
# - kern: array of weights for smoothing generated from the function kern.3d above
# Output: 3-dimensional array of the same size showing the result of applying 
# the smoother to the data
spatial.smooth3d <- function(image, kern)
{
  n1 <- dim(image)[1]; n2 <- dim(image)[2]; n3 <- dim(image)[3]
  kerndim <- dim(kern)[1]
  out <- .C('Rwrapper_spatial_smooth3d', as.double(image), 
            as.integer(n1), as.integer(n2), as.integer(n3), 
            as.double(kern), as.integer(kerndim), NAOK=T)[[1]]
  array(out, dim=c(n1, n2, n3))
}