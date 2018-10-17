Authors: Adrien Meynard and Bruno Torrésani  
Email: adrien.meynard@univ-amu.fr  
Created: 2017-12-19

This repository contains all the files necessary to perform the algorithm JEFAS (Joint Estimation of Frequency, Amplitude and Spectrum). The theoretical background relative to the algorithm can be found in the following paper:  
**A. Meynard and B. Torrésani, "Spectral Analysis for Nonstationary Audio", *IEEE/ACM Transactions on Audio,
Speech and Language Processing*, vol. 26, pp. 2371–2380, Dec. 2018, [available online](https://hal.archives-ouvertes.fr/hal-01670187).**

The implementation of the algorithm uses MATLAB functions from the Optimization Toolbox (e.g. `fmincon`). It is therefore necessary to have access to this toolbox in order to run JEFAS. A less efficient version of the algorithm that does not require the Optimization Toolbox is available upon demand.

Descriptions of the different folders:
- `cwt`: contains the functions to compute continuous wavelet transforms and inverse wavelet transforms using the wavelets given in the paper (equation (2)). A script `example_cwt.m` enabling the display of the sharp wavelet is also given (including Supplementary material, figs. 1 and 2). 
- `JEFAS_algo`: contains all the functions necessary to implement JEFAS algorithm (subfolder `estimation`), together with some functions enabling the analysis (subfolder `analysis`) of the results of the JEFAS algorithm (baseline estimations, CRLB, stationarization). It also contains functions enabling cross-synthesis.
- `signals`: contains some audio signals including those described in the article.
- `scriptsIEEE_TASLP`: contains the scripts detailed above, together with the published versions (subfolder `html`), and some corresponding results of JEFAS (time warping functions, amplitude modulation functions, and spectrum) in subfolder `results`.

Descriptions of the different MATLAB scripts in the folder scriptsIEEE_TASLP:
- `JEFAS_toysig.m`: estimation of the deformations and spectrum from a synthetic signal (see section IV.A of the article). This script generates Fig. 1, Fig. 2 and the values in Table I (of the paper).
- `JEFAS_dolphin.m`: estimation of the spectrum of the underlying stationary signal from a recording a a dolphin sound". This script generates Fig. 3.
- `JEFAS_dopplerf1.m`: JEFAS estimation from the recording of a Formula 1 (from a fixed location). JEFAS estimates the deformation due to the Doppler effect (see section IV.B). This script generates Fig. 4.
- `JEFAS_sing.m`: estimation of the deformation from the sound produced by a singing female voice. This script also allows you to listen to the sound obtained after an estimated "stationarization".
- `JEFAS_wind.m`: estimation of the deformation from a recording of wind howling (see section V in the Supplementary material). 
- `cross_synthesis.m`: synthesis of nonstationary sounds starting from a sound stationarized by JEFAS algorithm, and deforming it using time warping and amplitude modulation functions estimated from another signal (see section VI in the Supplementary material).

All the sounds and the published versions of the Matlab scripts are also [available online](http://meynard.perso.math.cnrs.fr/paperJEFAS/NonStationaryAudio.html).
