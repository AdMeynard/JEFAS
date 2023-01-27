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

import numpy as np

#%%

def stationarize(y,dgamma,a=None):
    """STATAMTW:  estimation of the underlying stationary signal
    usage:	z = stat(y,a,dgamma)
    
    Input:
      y: non-stationary signal
      a (optional): amplitude modulation function
      dgamma: derivative of the time warping function gamma
    Output:
      z : estimated stationary signal
  """
    
    if a==None:
        Dz = y
    else:
        Dz = y / a
        
    z = invwarp(Dz,dgamma)
    
    return z

def invwarp(y,dgamma):

    N = len(y)
    ts = np.linspace(0,1,N)
    
    gamma = np.cumsum(dgamma)
    gamma = (gamma - gamma[0])/(gamma[-1] - gamma[0]) # normalization such that gamma[-1] = 1
    
    invgamma = np.interp(ts,gamma,ts)
    invdgamma = np.interp(invgamma,ts,dgamma) # \gamma'(\gamma^{-1})
    
    x = (1/np.sqrt(invdgamma)) * np.interp(invgamma,ts,y)
    
    return x

#%%

def unstationarize(z,a,dgamma):
    """deformAMWP:  synthesis of a new nonstationary signal
    usage:	y = deformAMWP(z,a,dgamma)
    
    Input:
      z: stationary signal
      a: amplitude modulation function
      dgamma: derivative of the time warping function gamma
    Output:
      y : synthetic nonstationary signal
      """

    Td = len(dgamma)
    Tz = len(z)
    td = np.linspace(0,1,Td)
    tz = np.linspace(0,1,Tz)
    
    dgamma = np.interp(tz,td,dgamma) # matching sizes
    a = np.interp(tz,td,a) #matching sizes
    

    Dz = warp(z,dgamma)
    y = a * Dz

    return y

def warp(x,dgamma):
    
    N = len(x)
    ts = np.linspace(0,1,N)
    
    gamma = np.cumsum(dgamma)
    gamma = (gamma - gamma[0])/(gamma[-1] - gamma[0]) # normalization such that gamma[-1] = 1

    y = np.sqrt(dgamma) * np.interp(gamma,ts,x)
    
    return y
                               


