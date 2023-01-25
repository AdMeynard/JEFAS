Descriptions of the different MATLAB scripts in the folder scriptsIEEE_TASLP:
- `JEFAS_toysig.m`: estimation of the deformations and spectrum from a synthetic signal (see section IV.A of the article). This script generates Fig. 1, Fig. 2 and the values in Table I (of the paper).
- `JEFAS_dolphin.m`: estimation of the spectrum of the underlying stationary signal from a recording a a dolphin sound". This script generates Fig. 3.
- `JEFAS_dopplerf1.m`: JEFAS estimation from the recording of a Formula 1 (from a fixed location). JEFAS estimates the deformation due to the Doppler effect (see section IV.B). This script generates Fig. 4.
- `JEFAS_sing.m`: estimation of the deformation from the sound produced by a singing female voice. This script also allows you to listen to the sound obtained after an estimated "stationarization".
- `JEFAS_wind.m`: estimation of the deformation from a recording of wind howling (see section V in the Supplementary material). 
- `cross_synthesis.m`: synthesis of nonstationary sounds starting from a sound stationarized by JEFAS algorithm, and deforming it using time warping and amplitude modulation functions estimated from another signal (see section VI in the Supplementary material).
