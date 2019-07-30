# pointWav
Estimation of Smoothed Wavelet Spectra/Coherence in Continuous time from Multivariate Point-Processes

## Motivation

Observations from point-processes are often recorded in continuous time, binning these events can allow for analysis as a 
discrete time-series. However, often high-frequency dependency structure can be lost when performing such aggregations. To avoid
this, we have developed a temporal smoothing method for the continuous wavelet transform (and corresponding wavelet spectra).
The spectral estimators in this package thus operate on the continuous time event processes directly.

## How to Use

1. Simply add the folder and its contents to your MATLAB file path
2. Set up your data as a set of cells ```E{1}, E{2}, ... E{p}``` each with a vector of points
3. Generate spectral estimates ```S = tsWP( E, kappa, wavType ,'method','kernel');```
 * kappa =  smoothing window scale
 * wavType = Morlet, Mex (only these two are implemented currently)
