#%% Deformation estimations on a racing car sound
#
# Copyright (C) 2023 Adrien MEYNARD

# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Author: Adrien MEYNARD
# Email: adrien.meynard@ens-lyon.fr

#%%  Load signal
import matplotlib.pyplot as plt
from scipy.io import loadmat
import sys

sys.path.append(r'../JEFAS')

import JEFAS
import stationaryJEFAS
import cwtJEFAS

plt.close('all')

#%% Signal

tmp = loadmat('../signals/doppler_f1')
y = tmp['y'].flatten()
Fs = float(tmp['Fs'])

#%% JEFAS

Dt = 100 # hop size for the deformation estimation

wavparam = 20 # wavelet parameter for warping estimation

freqmin = 0.004
freqmax = 0.125
nf = 42

wavparamAM = 500

# Estimation
dgamma1, Sx1, evolcrit1 = JEFAS.JEFAS(y,Dt,wavparam,freqmin,freqmax,nf,disp=1,wavparamAM=wavparamAM)

z = stationaryJEFAS.stationarize(y, dgamma1)

plt.figure()
plt.subplot(121)
Wy,_,_ = cwtJEFAS.display_cwt(y, Fs, Fs*freqmin, Fs*freqmax,200,'sharp', 300)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.title('Original Signal')
plt.subplot(122)
Wz,t,s = cwtJEFAS.display_cwt(z, Fs, Fs*freqmin, Fs*freqmax,200,'sharp', 300)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.title('Stationarized Signal')
