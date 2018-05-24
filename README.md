Author: Adrien MEYNARD  
Email: adrien.meynard@univ-amu.fr  
Created: 2017-12-19

This repository contains all the files necessary to perform the algorithm JEFAS (Joint Estimation of Frequency, Amplitude and Spectrum). The theoretical background relative to the algorithm can be found in the paper "Spectral analysis for nonstationary audio" (A. Meynard and B. Torrésani), [available online](https://hal.archives-ouvertes.fr/hal-01670187).

The implementation of the algorithm uses MATLAB functions from the optimization toolbox. Thus, this toolbox is necessary to run the algorithm. A less efficient version of the algorithm that does not require this toolbox is available upon demand.

Description of the different MATLAB scripts:
- `JEFAS_toysig`: estimation of the deformation from the synthetic signal (see section IV.A). This script generates Fig. 1, Fig. 2 and the values in Table I.
- `JEFAS_dolphin`: estimation of the spectrum of the underlying stationary signal from a recording a a dolphin sound". This script generates Fig. 3.
- `JEFAS_dopperf1`: estimation of the deformation from the sound (see section IV.B). This script generates Fig. 4.
- `JEFAS_sing`: estimation of the deformation from the sound produced by a singing female voice. This script also allows you to listen to the sound obtained after an estimated "stationarization". 
- `cross_synthesis`: synthesis of nonstationary sounds starting from a sound stationarized by JEFAS algorithm, and deforming it using time-warping and amplitude modulation functions estimated from another signal.

All the sounds and the published versions of the Matlab scripts are [available online](http://meynard.perso.math.cnrs.fr/paperJEFAS/NonStationaryAudio.html).

Description of the different folders:
- html: this folder contains the published versions (html) of these three scripts.
- cwt: this folder contains the functions to compute continuous wavelet transforms and inverse wavelet transforms using the wavelets given in the paper (equation (2)). A script `example_cwt` enabling the display of the sharp wavelet is also given. 
- deform_estimation: this folder contains all the functions necessary to implement the alternate algorithm 1.
- analysis: this folder contains some functions enabling the analysis of the results of the algorithm JEFAS (baseline estimations, CRLB, stationarization). It also contains functions enabling cross-synthesis (!).
- signals: this folder contains some audio signals including those described in the article.
- results: this folder contains estimated time-warping functions, amplitude modulation fucntions, and spectrum for some examples.
