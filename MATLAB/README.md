

The algorithms require to implement the wavelet transforms of signals. The folder `cwt` contains the functions to compute continuous wavelet transforms and inverse wavelet transforms using the *sharp wavelet*, whose definition is given in [1] (see equation (2)). We also provide the script `example_cwt.m` enabling the display of the sharp wavelet and its Fourier transform (see Supplementary material of [1] for more details).

The folder `signals` contains some audio and synthetic signals used to implement the above-mentioned algorithms. Additionally, the folder `SynthSignals` contains functions to synthesize locally time-warped toy signals (both univariate and multivariate).

***WARNING:*** The implementation of these algorithms uses MATLAB functions from the Optimization Toolbox (e.g. `fmincon`). It is therefore necessary to have access to this toolbox in order to run JEFAS. A less efficient version of the algorithm that does not require the Optimization Toolbox is available upon demand.

# JEFAS
The following folders are related to the implementation of JEFAS:

* `JEFASalgo`: contains all the functions necessary to implement JEFAS (subfolder `estimation`), together with some functions enabling the analysis (subfolder `analysis`) of the results of JEFAS (comparison with baseline estimations, Cramér-Rao Lower Bound, stationarization). It also contains functions enabling cross-synthesis.
* `scriptsIEEE_TASLP`: contains the scripts detailed below, together with their published versions (subfolder `html`), and some corresponding results of JEFAS (estimated time warping functions, amplitude modulation functions, and spectra) in subfolder `results`. All the sounds and the published versions of the Matlab scripts are also [available online](http://meynard.perso.math.cnrs.fr/paperJEFAS/NonStationaryAudio.html).

# JEFAS-BSS
In addition to the folders described above, the following folders are related to the implementation of JEFAS-BSS:

* `JEFAS-BSS`: contains the functions to perform JEFAS-BSS, an algorithm for the spectral estimation of multivariate nonstationary signals following the model given in [1]. This problem can be equivalently seen as a doubly nonstationary Blind Source Separation (BSS) problem.
* `scriptsBSS`: contains the scripts that perform JEFAS-BSS on a synthetic mixture of both synthetic signals and real audio signals.

# JEFAS-S
In addition to the folders described above, the following folders are related to the implementation of JEFAS-S:

* `JEFAS-S`: contains the functions to perform JEFAS-S, a synthesis-based approach for the spectral estimation of signals following the model given in [1].
* `scriptsJEFAS-S`: contains the scripts that perform JEFAS-S on two examples.
