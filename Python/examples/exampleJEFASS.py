#%%% JEFAS-S on the narrowband signal

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import sys

sys.path.append(r'../JEFAS')

import JEFASS

plt.close('all')

tmp = loadmat('../signals/TwoSineWaves_fast.mat')
sigmay = tmp['sigmay']
dgamma = tmp['dgamma'].flatten()
Fs =float( tmp['Fs'])

tmp = loadmat('../signals/JEFASS_2sinesfast')
y = tmp['y'].flatten()

theta = np.log2(dgamma)

T = len(y)

#%% JEFAS-S

Nbscales = 18
scales = 2**np.linspace(-0.67,2.21,Nbscales)

wavparam = 20

TT = 512 #  slicing size for W

Nit = 50

dgammaEST, SxEST, W, nll = JEFASS.JEFASS(y,sigmay,scales,wavparam,TT,Nit=Nit,dgammaT=dgamma)


#%% Results

plt.figure()
plt.plot(nll)
plt.xlabel('Iteration') 
plt.ylabel('Negative log-likelihood of the signal')

t = np.arange(T)/Fs

plt.figure()
plt.plot(t,dgamma,'k--')
plt.plot(t,dgammaEST,'r')
plt.xlabel('Time (s)')
plt.ylabel(r"$\gamma'$")
plt.ylim(0.5,1.5)
