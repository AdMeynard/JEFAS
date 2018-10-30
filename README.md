Authors: Adrien Meynard and Bruno Torrésani  
Email: adrien.meynard@univ-amu.fr  
Created: 2017-12-19

This repository contains all the files necessary to perform the algorithm JEFAS (Joint Estimation of Frequency, Amplitude and Spectrum). The theoretical background relative to the algorithm can be found in the following paper:  
**[1] A. Meynard and B. Torrésani, "Spectral Analysis for Nonstationary Audio", *IEEE/ACM Transactions on Audio,
Speech and Language Processing*, vol. 26, pp. 2371–2380, Dec. 2018, [available online](https://hal.archives-ouvertes.fr/hal-01670187).**

The implementation of the algorithm uses MATLAB functions from the Optimization Toolbox (e.g. `fmincon`). It is therefore necessary to have access to this toolbox in order to run JEFAS. A less efficient version of the algorithm that does not require the Optimization Toolbox is available upon demand.

Descriptions of the different folders:
- `cwt`: contains the functions to compute continuous wavelet transforms and inverse wavelet transforms using the wavelets given in the paper (equation (2)). A script `example_cwt.m` enabling the display of the sharp wavelet is also given (including Supplementary material, figs. 1 and 2). 
- `JEFASalgo`: contains all the functions necessary to implement JEFAS algorithm (subfolder `estimation`), together with some functions enabling the analysis (subfolder `analysis`) of the results of the JEFAS algorithm (baseline estimations, CRLB, stationarization). It also contains functions enabling cross-synthesis.
- `signals`: contains some audio signals including those described in the article.
- `scriptsIEEE_TASLP`: contains the scripts detailed below, together with their published versions (subfolder `html`), and some corresponding results of JEFAS (estimated time warping functions, amplitude modulation functions, and spectra) in subfolder `results`. All the sounds and the published versions of the Matlab scripts are also [available online](http://meynard.perso.math.cnrs.fr/paperJEFAS/NonStationaryAudio.html).
- `JEFAS-BSS`: contains the functions to performs a Blind Source Separation (BSS) from a mixture of nonstationary signals following the model given in [1].
- `scriptsIEEE_TASLP`: contains the scripts that perform JEFAS-BSS on a synthetic mixture.
