# % Deformation estimations on a singing female voice recording
# Copyright (C) 2022 Adrien MEYNARD

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
# Email: adrien.meynard@univ-amu.fr

#%%  Load signal
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
import sys

sys.path.append(r'../JEFAS')

import cwtJEFAS

plt.close('all')

tmp = loadmat('../../signals/sing')
y = tmp['y']
Fs = float(tmp['Fs'])

#%%
N = len(y)
t = np.arange(N)/Fs

wav_typ = 'sharp'
wav_param = 100

plt.figure()
freqtick = [200, 450, 1000]
Wy,_,scales = cwtJEFAS.display_cwt(y,Fs,200,2000,150,wav_typ,wav_param,freqtick=freqtick)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
plt.title('CWT')

yr = cwtJEFAS.icwt(Wy,scales,wav_typ,wav_param)

plt.figure()
plt.plot(t,y,label='Original signal')
plt.plot(t,yr,label='Synthesized signal')
plt.xlim(0,0.2)
plt.xlabel('Time (s)')
plt.ylabel('Signals')
plt.legend()
plt.title('Signal synthesis from CWT')