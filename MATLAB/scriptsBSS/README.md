***WARNING:*** Make sure to run the file `download_data.m` once before launching the script `JEFASBSS_SynthMixs.m` or `JEFASBSS_SynthMixsOverDeterm.m`. This allows you to download the required data to run these scripts.

Descriptions of the different MATLAB scripts in the folder scriptsBSS:
- `JEFASBSS_SynthMix.m`: performs JEFAS-BSS on one realization of a determined synthetic mixture of nonstationary synthetic signals. The performance of JEFAS-BSS is compared with baseline BSS algorithms.
- `JEFASBSS_SynthMixs.m`: performs JEFAS-BSS on 20 realizations of the same synthetic mixture. Two performances indexes are evaluated and compared with baseline BSS algorithms.
- `JEFASBSS_SynthMixOverDeterm.m`: performs JEFAS-BSS on one realization of an overdetermined synthetic mixture nonstationary synthetic signals i.e. the number of observations is larger than the number of sources. The performance of JEFAS-BSS are compared with baseline BSS algorithms.
- `JEFASBSS_SynthMixsOverDeterm.m`: performs JEFAS-BSS on 20 realizations of the same overdetermined synthetic mixture. Two performance indexes are evaluated and compared with baseline BSS algorithms.
- `JEFASBSS_SynthMixOverDetermIllCond.m`: performs JEFAS-BSS on one realization of an overdetermined synthetic mixture nonstationary synthetic signals in a situation where one the mixing submatrices is ill-conditioned. The results provided by this script allow us to highlight the fact that JEFAS-BSS outperforms *Single JEFAS-BSS* in that case.
- `JEFASBSS_sounds.m`: performs JEFAS-BSS on a synthetic mixture of real nonstationary audio sources.

Descriptions of the different subforders:
- `CompSpeedMat`: contains MATLAB scripts enabling the evaluation of the BSS performance in function of the speed of variation of the mixing matrix coefficients.
- `results`: contains some BSS results.
