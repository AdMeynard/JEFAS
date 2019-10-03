***WARNING:*** Make sure to run the file `download_data.m` once before launching the script `JEFASBSS_SynthMixs.m` or `JEFASBSS_SynthMixsOverDeterm.m`. This allows you to download the required data to run these scripts.

Descriptions of the different MATLAB scripts in the folder scriptsBSS:
- `JEFASBSS_SynthMix.m`: performs JEFAS-BSS on one realization of a determined synthetic mixture of nonstationary synthetic signals. The performances of JEFAS-BSS are compared with baseline BSS algorithms.
- `JEFASBSS_SynthMixs.m`: performs JEFAS-BSS on 20 realizations of the same synthetic mixture. Two performances indexes are evaluated and compared with baseline BSS algorithms.
- `JEFASBSS_SynthMixOverDeterm.m`: performs JEFAS-BSS on one realization of an overdetermined synthetic mixture nonstationary synthetic signals i.e. the number of observations is larger than the number of sources. The performances of JEFAS-BSS are compared with baseline BSS algorithms.
- `JEFASBSS_SynthMixsOverDeterm.m`: performs JEFAS-BSS on 20 realizations of the same overdetermined synthetic mixture. Two performances indexes are evaluated and compared with baseline BSS algorithms.
- `JEFASBSS_sounds.m`: performs JEFAS-BSS on a synthetic mixture of real nonstationary audio sources.
