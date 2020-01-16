# Contents
This repository contains all the files necessary to perform the following algorithms:

* **JEFAS**: Joint Estimation of Frequency, Amplitude and Spectrum
* **JEFAS-BSS**: JEFAS-Blind Source Separation
* **JEFAS-S**: JEFAS-Synthesis

***WARNING:*** The implementation of these algorithms uses MATLAB functions from the Optimization Toolbox (e.g. `fmincon`). It is therefore necessary to have access to this toolbox in order to run JEFAS. A less efficient version of the algorithm that does not require the Optimization Toolbox is available upon demand.

# JEFAS
The theoretical background relative to the algorithm can be found in the following paper:  
**[1] A. Meynard and B. Torrésani, "Spectral Analysis for Nonstationary Audio", *IEEE/ACM Transactions on Audio,
Speech and Language Processing*, vol. 26, pp. 2371–2380, Dec. 2018, [available online](https://hal.archives-ouvertes.fr/hal-01670187).**

The following folders are related to the implementation of JEFAS:

* `cwt`: contains the functions to compute continuous wavelet transforms and inverse wavelet transforms using the wavelets given in the paper (equation (2)). A script `example_cwt.m` enabling the display of the sharp wavelet is also given (including Supplementary material, figs. 1 and 2). 
* `JEFASalgo`: contains all the functions necessary to implement JEFAS (subfolder `estimation`), together with some functions enabling the analysis (subfolder `analysis`) of the results of JEFAS (comparison with baseline estimations, Cramér-Rao Lower Bound, stationarization). It also contains functions enabling cross-synthesis.
* `signals`: contains some audio signals including those described in the article.
* `scriptsIEEE_TASLP`: contains the scripts detailed below, together with their published versions (subfolder `html`), and some corresponding results of JEFAS (estimated time warping functions, amplitude modulation functions, and spectra) in subfolder `results`. All the sounds and the published versions of the Matlab scripts are also [available online](http://meynard.perso.math.cnrs.fr/paperJEFAS/NonStationaryAudio.html).

# JEFAS-BSS
The theoretical background relative to the algorithm can be found in the following paper:  
**[2] A. Meynard, "Spectral Estimation for Multivariate Locally Time-Warped Signals", submitted**

In addition to the folders described above, the following folders are related to the implementation of JEFAS-BSS:

* `JEFAS-BSS`: contains the functions to perform JEFAS-BSS, an algorithm for the spectral estimation of multivariate nonstationary signals following the model given in [1]. This problem can be equivalently seen as a doubly nonstationary Blind Source Separation (BSS) problem.
* `scriptsBSS`: contains the scripts that perform JEFAS-BSS on a synthetic mixture of both synthetic signals and real audio signals.

# JEFAS-S
The theoretical background relative to the algorithm can be found in Chapter 4 of the following thesis:  
**[3] A. Meynard : "Stationnarités brisées : approches à l'analyse et à la synthèse", PhD thesis, Aix-Marseille Université, Oct. 2019**

In addition to the folders described above, the following folders are related to the implementation of JEFAS-S:

* `JEFAS-S`: contains the functions to perform JEFAS-S, a synthesis-based approach for the spectral estimation of signals following the model given in [1].
* `scriptsJEFAS-S`: contains the scripts that perform JEFAS-S on two examples.

# Authors

Authors: Adrien Meynard and Bruno Torrésani  
Contact email: adrien.meynard@duke.edu  
Created: 2017-12-19

