Author: Adrien MEYNARD  
Email: adrien.meynard@univ-amu.fr  
Created: 2017-12-19

This repository contains all the files necessary to perform relative to the publication "Spectral analysis for nonstationary audio" (A. Meynard and B. Torrésani), [available here](https://hal.archives-ouvertes.fr/hal-01670187)


The implementation of the algorithm uses MATLAB functions from the optimization toolbox. Thus, this toolbox is necessary to run the algorithm. A less efficient version of the algorithm that does not require this toolbox is available upon demand.


Description of the different MATLAB scripts:
- `estimWPAM_toysig`: estimation of the deformation from the synthetic signal (see section IV.A). This script generates Figure 1, Figure 2 and the values in Table I.
- `estimWPAM_dopperf1`: estimation of the deformation from the sound (see section IV.B). This script generates Figure 3.
- `estimWPAM_sing`: estimation of the deformation from the sound produced by a singing female voice. This script also allows you to listen to the sound obtained after an estimated "stationarization".


Description of the different folders:
- html: this folder contains the published versions (html) of these three scripts.
- cwt: this folder contains the functions to compute continuous wavelet transforms and inverse wavelet transforms with the wavelets given in the paper (equation (2)).
- deform_estimation: this folder contains all the functions necessary to implement the alternate algorithm 1.
- analysis: this folder contains some functions enabling the analyze the results of the algorithm 1 (baseline estimations, CRLB, stationarization).
- signals: this folder contains some audio signals including those described in the article.
