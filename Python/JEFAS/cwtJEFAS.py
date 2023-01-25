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

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

#%% CWT

def cwt(sig,scales,wavtyp,wavpar):
    """cwt:	Continuous wavelet transform with analytic derivative of Gaussian
          wavelet or the sharp wavelet
    usage:	W = cwt(sig,scales,wavtyp,wavpar)
    
    Input:
      sig: vector containing the input signal
      scales : vector of scales
      wavtyp: wavelet type:
        wavtyp='sharp': The sharp wavelet
        wavtyp='dgauss': analytic derivative of gaussian
      wavpar: parameter depending on the wavelet type
        if wavtyp='sharp': wavpar = -ln(epsilon)>0 where epsilon=value of \hat{\psi}(Fs/2)
        if wavtyp='dgauss': wavpar = number of vanishing moments
    
    Output:
      W: wavelet transform coefficients matrix (complex valued)
      """

    siglength = len(sig)
    Nsc = len(scales)
    
    fff = np.arange(siglength)*2*np.pi/siglength
    tmp = np.outer(scales, fff)
    
    scales = np.reshape(scales,(Nsc,1)) # column vector
    
    sig = np.reshape(sig,(1,siglength)) # row vector
    fsig = np.fft.fft(sig)
    
    
    # Generate wavelets in the Fourier domain
    if wavtyp=='sharp':
        eps = np.finfo(np.float64).eps # numpy precision
        fpsi = np.exp(-2*wavpar*((np.pi/(2*tmp+eps)) + (2*tmp/np.pi) - 2))
    
    elif wavtyp=='dgauss':
        Cst = 4*wavpar/(np.pi**2) # such that argmax fpsi = pi/2
        K = (2/np.pi)**wavpar*np.exp(wavpar/2) #such that fpsi(pi/2) = 1
        fpsi = K * tmp**wavpar * np.exp(-Cst*tmp**2/2)
        
    else:
        raise NameError('Unexpected wavelet type. CWT not computed')    
    
    U = np.sqrt(scales) * fpsi
    fTrans = U * fsig
    W = np.fft.ifft(fTrans)
    
    return W


#%% Display function

def display_cwt(sig,Fs,fmin,fmax,Ms,wavtyp,wavpar,freqtick=None,disp=True):
    """display_cwt:	Display the continuous wavelet transform with analytic derivative of Gaussian
          wavelet or the sharp wavelet
    usage:	W = display_cwt(sig,Fs,fmin,fmax,Ms,wavtyp,wavpar,freqtick,disp=True)
      /!\ If you do not want to display the CWT, set disp=False
    
    Input:
      sig: vector containing the input signal
      Fs: sampling frequency
      fmin: minimum frequency to display
      fmax: maximum frequency to display
      Ms: length of the scales vector (= number of rows of the CWT) 
      wavtyp: wavelet type:
        wavtyp='sharp': The sharp wavelet
        wavtyp='dgauss': analytic derivative of gaussian
      wavpar: parameter depending on the wavelet type
        if wavtyp='sharp': wavpar = -ln(epsilon)>0 where epsilon=value of \hat{\psi}(Fs/2)
        if wavtyp='dgauss': wavpar = number of vanishing moments
      freqtick (optional): frequency to tick on the y-axis
    
    Output:
      W: wavelet transform coefficients matrix (complex valued)
      """

    N = len(sig)
    t = np.arange(N)/Fs
    
    # Scales
    xi0 = Fs/4 # wavelet central frequency
    smin = np.log2( xi0/fmax )
    smax = np.log2( xi0/fmin )
    scales = 2**(np.linspace(smin,smax,Ms))
    
    # Wavelet transform:
    W = cwt(sig,scales,wavtyp,wavpar)
    
    
    # Display:
    if not disp :
        return W,t,scales
    else:
        if freqtick==None:
            freqtick = np.floor(fmax*(fmin/fmax)**(np.linspace(0,1,5)))
        else:
            freqtick = np.array(freqtick)
            
        sdisp = np.log2(xi0/freqtick) # corresponding log-scales
        
        lims = [t[0] , t[-1], smax , smin]
    
        plt.imshow(np.log1p(np.abs(W))/0.1,aspect='auto',extent=lims)
        plt.yticks(sdisp,freqtick)
        plt.show()
        return W,t,scales


#%% Inverse CWT
def icwt(W,scales,wavtyp,wavpar):
    """ICWT:	Continuous inverse wavelet transform with analytic derivative of
          Gaussian or the sharp wavelet
    usage:	y = icwt(W,scales,wavtyp,wavpar)
    
    Input:
      W: wavelet transform coefficients matrix
      scales : vector of scales
      wav_typ: wavelet type:
        wavtyp='sharp': The sharp wavelet
        wavtyp='dgauss': analytic derivative of gaussian
      wavpar: parameter depending on the wavelet type
        if wavtyp='sharp': wav_par = -ln(epsilon)>0 where epsilon=value of \hat{\psi}(Fs/2)
        if wavtyp='dgauss': wav_par = value of \hat{\psi}(Fs/2)
    
    Output:
      y: reconstructed signal
      """


    Nsc, siglength = np.shape(W)
    fff = np.arange(siglength)*2*np.pi/siglength
    
    scales = np.reshape(scales,(Nsc,1))
    ds = np.log(scales[1]) - np.log(scales[0])
    
    tmp = np.outer(scales, fff)
    
    # Generate wavelets in the Fourier domain
    if wavtyp=='sharp':
        eps = np.finfo(np.float64).eps # numpy precision
        fpsi = np.exp(-2*wavpar*((np.pi/(2*tmp+eps)) + (2*tmp/np.pi) - 2))
        Cpsi = 0.88/wavpar**0.5 # empirical
    
    elif wavtyp=='dgauss':
            Cst = 4*wavpar/(np.pi**2) # such that argmax fpsi = pi/2
            K = (2/np.pi)**wavpar*np.exp(wavpar/2) # such that fpsi(pi/2) = 1
            fpsi = K * tmp**wavpar * np.exp(-Cst*tmp**2/2)
            Cpsi = (K**2/(2*Cst**wavpar))*gamma(wavpar) # theoretical
        
    else:
       raise NameError('Unexpected wavelet type. ICWT not computed')    
    
    fy = fpsi * np.fft.fft(W)
    fy = 2*np.fft.ifft(fy).real
    fy = 1/scales**0.5 * fy
    y = (1/Cpsi) * np.sum(fy,0) * ds
    
    return y


    