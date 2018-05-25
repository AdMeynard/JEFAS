Authors: Adrien Meynard and Bruno Torrésani  
Email: adrien.meynard@univ-amu.fr  
Created: 2017-12-19

This repository contains all the files necessary to perform the algorithm JEFAS (Joint Estimation of Frequency, Amplitude and Spectrum). The theoretical background relative to the algorithm can be found in the paper "Spectral analysis for nonstationary audio" (A. Meynard and B. Torrésani), [available online](https://hal.archives-ouvertes.fr/hal-01670187).

The implementation of the algorithm uses MATLAB functions from the optimization toolbox. Thus, this toolbox is necessary to run the algorithm. A less efficient version of the algorithm that does not require this toolbox is available upon demand.

Description of the different MATLAB scripts:
- `JEFAS_toysig.m`: estimation of the deformation from the synthetic signal (see section IV.A). This script generates Fig. 1, Fig. 2 and the values in Table I.
- `JEFAS_dolphin.m`: estimation of the spectrum of the underlying stationary signal from a recording a a dolphin sound". This script generates Fig. 3.
- `JEFAS_dopperf1.m`: estimation of the deformation from the sound (see section IV.B). This script generates Fig. 4.
- `JEFAS_sing.m`: estimation of the deformation from the sound produced by a singing female voice. This script also allows you to listen to the sound obtained after an estimated "stationarization". 
- `cross_synthesis.m`: synthesis of nonstationary sounds starting from a sound stationarized by JEFAS algorithm, and deforming it using time-warping and amplitude modulation functions estimated from another signal.

All the sounds and the published versions of the Matlab scripts are [available online](http://meynard.perso.math.cnrs.fr/paperJEFAS/NonStationaryAudio.html).

Description of the different folders:
- html: contains the published versions (html) of these scripts.
- cwt: contains the functions to compute continuous wavelet transforms and inverse wavelet transforms using the wavelets given in the paper (equation (2)). A script `example_cwt.m` enabling the display of the sharp wavelet is also given. 
- deform_estimation: contains all the functions necessary to implement JEFAS algorithm.
- analysis: contains some functions enabling the analysis of the results of the JEFAS algorithm (baseline estimations, CRLB, stationarization). It also contains functions enabling cross-synthesis (!).
- signals: contains some audio signals including those described in the article.
- results: contains estimated time-warping functions, amplitude modulation fucntions, and spectrum for some examples.
