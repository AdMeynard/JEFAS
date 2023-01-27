#%% Comparison between time-scale representations of the narrowband signal

#%% Load signal and JEFAS-S results

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import sys

sys.path.append(r'../JEFAS')

import cwtJEFAS
import JEFASS

plt.close('all')

tmp = loadmat('../signals/TwoSineWaves_fast.mat')
sigmay = tmp['sigmay']

tmp = loadmat('../signals/JEFASS_2sinesfast')
y = tmp['y']
SxEST = tmp['SxEST'].flatten()
dgammaEST = tmp['dgammaEST'].flatten()

T = len(y)
Fs = T

#%% Show CWT

Ms = 200
fmin = np.floor(0.048*T)
fmax = np.floor(0.40*T)

wavtyp = 'sharp' # wavelet type
wavparam = 50
freqtick = [50, 100, 200, 300, 400]

plt.figure()

Wy,t,scales = cwtJEFAS.display_cwt(y,Fs,fmin,fmax,Ms,wavtyp,wavparam,freqtick=freqtick)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')

#%% Show ordinary Synthesis Time-Scale Representation

prior = 'sharp'
W, MMSigmay = JEFASS.TimeScaleRep(y,sigmay,prior,scales,wavparam,SxEST,dgammaEST)


freqtick = np.array(freqtick)
xi0 = Fs/4 # wavelet central frequency
sdisp = np.log2(xi0/freqtick) # corresponding log-scales

smin = np.log2( xi0/fmax )
smax = np.log2( xi0/fmin )
lims = [t[0] , t[-1], smax , smin]

plt.figure()
plt.imshow(np.log1p(np.abs(W)/0.1),aspect='auto',extent=lims)
plt.yticks(sdisp,freqtick)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')

