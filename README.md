# complex_fmri
This repo contains R and C code accompanying the Annals of Applied Statistics paper 
["Complex-valued time series modeling for improved activation detection in fMRI studies"](https://projecteuclid.org/euclid.aoas/1536652961).

The R functions are intended to be the user interface so I have written comments about them.  The C functions are called by these R functions and represent the "inner workings" of the code so I have not commented these (though I am happy to include more upon request).


Here are summaries of what the functions do in each of the R files:

* `processing-R_fcns.R` Performs the detrending of the magnitude/phase time series components via the complex-valued running line as described in Section 3.1; performs spatial smoothing of the fMRI images as described in Section 3.2.
* `order_det_fcns.R` Detects the autoregressive (AR) order of the complex-valued time series as described in Section 3.5.
* `complex_fcns.R` For a complex-valued voxel time series, estimates the parameters and calculates the likelihood ratio test (LRT) statistic for activation.

I have only included the code related to the models for complex-valued fMRI time series here.  I have excluded the code related to the models for the magnitude-only time series, though I can include this upon request.
