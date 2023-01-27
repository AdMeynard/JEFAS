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
from scipy.optimize import fmin
from scipy.signal import welch

import stationaryJEFAS

#%%

def JEFASS(y,sigmay,scales,wavparam,TT,wavtyp='sharp',Dt=1,alpha=7,Nit=100,errtol=1e-4,disp=1,**kargs):
    """EMWARPING	JEFAS-S joint estimation of the time warping function, spectrum, and time-scale representation
    
    Parameters
    ----------
    y : array_like
        Signal to analyze
    sigmay: float
        Noise variance
    scales : array_like
        Scales on which the time-scales representation is estimated        
    wavparam : float 
        Wavelet parameter for the time warping estimation (cf. cwtJEFAS)
    TT: int
        slicing size for the time-scale representation estimation
    wavtyp : string, optional
        Wavelet type (cf. cwtJEFAS), default is wavtyp='sharp'    
    Dt : int, optional
        Subsampling step for estimations (>=1), default is 1
    alpha: float
        Width for Welch estimation        
    Nit : int, optional
        Maximum number of iterations of the alternate algorithm (default: Nit = 10)
    errtol : float, optional
        Relative error in TW, AM and spectrum estimates between iterations that is acceptable for convergence (default is errtol = 0.01%)
    disp : int, optional
        If disp = 0, disable relative update evolution (default: disp = 1)
        
    Other Parameters
    ----------
    TWinit : array_like, optional
        Initial estimation of gamma'(t), default is TWinit[t]=1 for all t
    TS : int, optional
        Length of the spectrum estimator
    Delta : float, optional
        overlap for the slices in the TS representaion estimation
    gammaT : array_like, optionnal
        Ground-thruth TW function
      
    Returns
    ----------
    dgammaEST : ndarray
        Estimation of the time warping function
    SxEST : ndarray
        Estimation of the spectrum of the underlying stationary signal x
    W : ndarray
        Estimation of the time-scale representation
    nllV : ndarray
        Evolution of the signal negative log-likekihood through iterations
    """
    
    T = len(y)
    
    if 'TS' in kargs:
        TS = kargs['TS']
    else:
        TS=T
    
    if 'Delta' in kargs:
        Delta = round(kargs['Delta'])
    else:
        Delta = round(0.75*TT)
        
    prior = 'wavelet'
    
    nwin = np.floor(T/alpha)
    
    Mpsi = bas_calcCov(scales,wavtyp,wavparam,T)
    MatPsi = np.fft.ifft(Mpsi).T # Synthesis coefficients
    
    Mpsi = bas_calcCov(scales,wavtyp,wavparam,TS)
    
    if 'TWinit' in kargs:
        thetaold = np.log2(kargs['TWinit'])
    else:
        thetaold = np.zeros(T)
        
    if 'dgammaT' in kargs:
        dgammaT = kargs['dgammaT']
        thetaT = np.log2(dgammaT)
        errINIT = np.sum( np.abs( thetaold - thetaT)**2 )
        if disp==1:
            print('Initialization')
            print(f"Quadratic error: {errINIT:.3f}")
            print('')
    
    nbit = 1
    nll = 1/np.finfo(np.float64).eps
    errEM = np.inf
    errTW = np.inf
    
    nllV = []
    # EM alternate estimation
    while ( (nbit<=Nit) and (errTW>errtol) ):
        
        # Step E: Adapted transform estimation
        x = stationaryJEFAS.stationarize(y,2**thetaold) # no AM estimation here
        Sxest = welch(x,nperseg=nwin,nfft=TS,return_onesided=False)[1]
        W, MMSigmay = bas_TimeScaleRep(y,sigmay,prior,TT,Delta,Mpsi,Sxest,thetaold,MatPsi,scales=scales)
        Sigmay = buildSigmay(MMSigmay,T,TT,Delta) # full signal covariance matrix
        iSigmay = np.linalg.inv(Sigmay) # INVERSION !! DANGER !!
        
        # Step M: Time-warping estimation
        thetaEM = []
        for t in range(0,T,Dt):
            U = W[:,t]
            theta0 = thetaold[t]
    
            theta = fmin(Qem,theta0,args=(theta0, t, U, Mpsi, Sxest, MatPsi, iSigmay),disp=0) # Minimized function
            thetaEM = np.append(thetaEM,theta)
            
        thetanew = np.interp(np.arange(T),np.arange(0,T,Dt),thetaEM,left=thetaEM[0],right=thetaEM[-1]) # interpolation on all the samples
        
        
        nll = negloglikelihoodsig(y,iSigmay) # negative log-likelihood (decreasing)
        nllV = np.append(nllV, nll)
        
        if np.any(thetaold):
            errTW = np.sum((thetanew - thetaold)**2)/np.sum(thetaold**2)
        else:
            errTW = np.inf

        if disp==1:
            if 'dgammaT' in kargs:
                errEM = np.sum( np.abs( thetanew - thetaT )**2 )
                print(f'Iteration {nbit}')
                print(f'Quadratic Error: {errEM:.3f}')
            else:
                print(f'Iteration {nbit}')
                print(f'Negative loglikelihood: {nll:.3f}')
            
            print(f"Relative update TW: {100*errTW:.2f} %")
            print('')
            
        
        nbit = nbit + 1
        thetaold = thetanew
    
    dgammaEST = 2**thetaold 
    x = stationaryJEFAS.stationarize(y,dgammaEST)
    SxEST = welch(x,nperseg=nwin,nfft=TS,return_onesided=False)[1] # spectrum
    W = bas_TimeScaleRep(y,sigmay,prior,TT,Delta,Mpsi,Sxest,thetaold,MatPsi,scales=scales)
    
    Sigmay = buildSigmay(MMSigmay,T,TT,Delta) # full signal covariance matrix
    iSigmay = np.linalg.inv(Sigmay) # INVERSION !! DANGER !!
    nll = negloglikelihoodsig(y,iSigmay) # negative log-likelihood (decreasing)
    nllV = np.append(nllV, nll)
    nllV = nllV[1:]
    
    return dgammaEST,SxEST, W, nllV

#%%

def TimeScaleRep(y,sigmay,prior,scales,wavparam,Sx,dgamma,TT=None,Delta=None):
    
    T = len(y)

    if TT is None:
        TT = int(np.minimum(1024,np.round(T/2)))
    
    if Delta is None:
        Delta = int(0.75*TT)


    Mpsi = bas_calcCov(scales,'sharp',wavparam,T)
    MatPsi = np.fft.ifft(Mpsi).T # Synthesis coefficients

    if prior=='sharp':
        sigma_sharp = np.log2(scales[1]/scales[0])
        
        s = np.log2(scales)
        Ns = len(s)
        S = np.tile(s,(Ns,1))
        Mpsi = np.exp(-(S-S.T)**2/(2*sigma_sharp**2))
        
    theta = np.log2(dgamma)
        
    W,MMSigmay = bas_TimeScaleRep(y,sigmay,prior,TT,Delta,Mpsi,Sx,theta,MatPsi,scales=scales)
        
    return W,MMSigmay

def bas_TimeScaleRep(y,sigmay,prior,TT,Delta,Mpsi,Sx,theta,MatPsi,scales=None):
    """ bas_TimeScaleRep	determine the adapted transform of a signal given the warping function and underlying spectrum
     usage:	W, MMSigmay =  bas_TimeScaleRep(y,sigmay,prior,TT,Delta,Mpsi,Sx,theta,MatPsi,scales) 
    
     Input:
       y: signal
       sigmay: noise variance
       prior: shape of the a priori covariance matrix ('wavelet' or 'sharp')
       TT: slicing size for W
       Delta: overlap
       Mpsi: base for covariance computation
       Sx: underlying spectrum
       theta: WP parameters
       MatPsi: characterize the transform shape (wavelet)
       scales: vector of scales (only when prior!='wavelet' because scales is included in Mpsi)
     
     Output:
       W : adapted transform
       MMSigmay : basis for signal covariance matrix estimation
    """
    
    if prior=='wavelet' and scales is None:
        raise NameError('scales vector must be speficified')
    
    T = len(y)
    
    TS = len(Sx)
    omega = np.arange(TS)*2*np.pi/TS
    
    delta = int(np.floor( (TT-Delta)/2 ))
    K = int(np.ceil( (T-TT)/Delta ) + 1) # number of blocks
    
    MN = list(range(TT))
    for dec in range(TT):
        Mdec = np.roll(MatPsi,dec,axis=0)
        MN[dec] = Mdec[0:TT,:]
    
    nMax = -1
    MMSigmay = list(range(K))
    C = list(range(T))
    CMn = list(range(TT))
    W = np.zeros((len(scales),T),dtype=complex)
    
    for k  in range(K): # k-th block
        nn = np.arange( k*Delta, k*Delta+TT )
        nn = nn[(nn>=0)&(nn<T)] # time of the k-th block
        NN = len(nn)
        ycourt = y[nn] # we divide the signal into segments of size TT
        
        if np.isreal(MatPsi).all():
            Sigmay = sigmay**2*np.eye(NN)
        else:
            Sigmay = 2*sigmay**2*np.eye(NN)
        
        for n in nn:
            Mn = MN[n-min(nn)]
            if NN != TT:
                Mn = Mn[0:NN,:]
            if n > nMax: # avoid redundancies 
                if prior=='wavelet':
                    C[n] = calcCov_synth(prior,Mpsi,Sx,omega,theta[n],T=TS)
                elif prior=='sharp':
                    C[n] = calcCov_synth(prior,Mpsi,Sx,omega,theta[n],scales=scales)
                    
            CMn[n-min(nn)] = C[n] @ np.conj(Mn.T)
            Sigmay = Sigmay + Mn @ CMn[n-min(nn)]

        Sigmay = Sigmay.real
        MMSigmay[k] = Sigmay # we store successive Sigmay matrices
        
        # Adapted representation computation
        invSy = np.linalg.inv(Sigmay) @ ycourt
        if k == 0: # first block
            for n in nn:
                W[:,n] = (CMn[n-min(nn)] @ invSy).flatten()
        elif k<(K-1):
            nnC = nn[ delta:(TT-delta) ]
            for n in nnC:
                W[:,n] = (CMn[n-min(nn)] @ invSy).flatten()
        else:
            nnC = nn[ delta: ]
            for n in nnC:
                W[:,n] = (CMn[n-min(nn)] @ invSy).flatten()
        
        nMax = nn[-1]
    
    return W, MMSigmay

#%%

def calcCov_synth(prior,Mpsi,S,omega,theta,T=None,scales=None):
    """CALC_COV_SYNTH Computation the covariance matrix
     usage:	C = calc_cov_synth(prior,Mpsi,S,omega,theta,T,scales)
    
     Input:
       prior: shape of the a priori covariance matrix ('wavelet' or 'sharp')
       Mpsi: matrix output of the function BAS_CALC_COV
       S: spectrum of X (row vector)
       omega: vector of frequencies
       theta: Time warping parameter
       T: signal length (if prior=="wavelet")
       scales: vector of scales (if prior=="sharp")
    
     Output:
       C: covariance matrix
    """
    
    if prior=='wavelet':
        if T is None:
            raise ValueError('T must be specified')
        else:
            Stheta = np.interp(2**(-theta)*omega,omega,S,left=0,right=0)
            U = np.sqrt(Stheta) * Mpsi
            C = (U @ np.conj(U.T)) / T
        
    elif prior=='sharp':
        if scales is None:
            raise ValueError('scales must be specified')
        else:
            omega0 = np.pi/2 # the wavelet mode is at Fs/4
            vtheta = np.sqrt( np.interp(2**(-theta)*omega0/scales,omega,S,left=0,right=0) ) 
            C = np.outer(np.conj(vtheta), vtheta) * Mpsi    
         
    return C

#%%

def bas_calcCov(scales,wavtyp,wavparam,T):

    omega = np.arange(T)*2*np.pi/T
    
    M_tmp = np.outer(scales,omega)
    scales = np.reshape(scales,(len(scales),1))
    
    if wavtyp=='sharp':
        eps = np.finfo(np.float64).eps # numpy precision
        par = -2*wavparam
        Mpsi = np.exp(par*((np.pi/(2*M_tmp+eps)) + (2*M_tmp/np.pi) - 2))
    elif wavtyp=='dgauss':
        Cst = 4*wavparam/(np.pi**2)
        K = (2/np.pi)**wavparam*np.exp(wavparam/2)
        Mpsi = K * M_tmp**wavparam * np.exp(-Cst*M_tmp**2/2)
    
    Mpsi = np.sqrt(scales) * Mpsi
    
    return Mpsi

#%%

def buildSigmay(MMSigmay,T,TT,Delta):
    """BUILDSIGMAY Concatenate short covariance matrices of subsignals to
    create the full covariance matrix of the signal
     usage:	Sigmay = buildSigmay(MMSigmay,T,TT,Delta)
    
     Input:
       MMSigmay: basis for signal covariance matrix estimation
       T: signal length
       TT: slicing size for W
       Delta: overlap
     
     Output:
       Sigmay: signal covariance matrix
    """
    
    K = len(MMSigmay)
    
    Sigmay = np.zeros((T,T))
    for k in range(K): # k-th time block
        nn = np.arange( k*Delta , k*Delta+TT )
        nn = nn[(nn>=0)&(nn<T)] # instants of the k-th block
        
        nn0 = nn[0]
        nnf = nn[-1]+1
        Sigmay[nn0:nnf,nn0:nnf] = MMSigmay[k]
    
    return Sigmay

#%%

def negloglikelihoodsig(y,iSigmay):
    """NEGLOGLIKELIHOODSIG negative log likelihood of the signal
     usage:	nllk = negloglikelihoodsig(y,iSigmay)
    
     Input:
       y: signal
       iSigmay: inverse signal covariance matrix
    
     Output:
       nllk: value of the negative log likelihood
    """
    
    nllk = np.dot(y, iSigmay @ y) - np.linalg.slogdet(iSigmay)[1] # Should decrease through iterations in EM algorithm
    
    return nllk

#%%

def Qem(theta, thetaprec, t, U, Mpsi, S, MatPsi, iSigmay):
    """Qem:  evaluation of the function Q to minimize during EM algorithm 
     usage:	[g,dg] = Qem(theta, thetaprec, t, U, Mpsi, S, MatPsi, iSigmay)
    
     Input:
       theta: value of the parameter
       thetaprec: value of the parameter given by the previous EM iteration
       t: time
       U: column of the tim-scale representation at time t
       Mpsi: first matrix output of the function BAS_CALC_DCOV
       S: spectrum
       MatPsi: matrix of "wavelets" for synthesis 
       iSigmay: inverse of the covariance matrix of the signal for synthesis
     Output:
       g : value of Q(theta, thetaprec)
    """
    
    T = len(S)
    omega = np.arange(T)*2*np.pi/T
    
    Cold = calcCov_synth('wavelet',Mpsi,S,omega,thetaprec,T=T) # covariance old
    C = calcCov_synth('wavelet',Mpsi,S,omega,theta,T=T) # covariance new

    Mn = np.roll(MatPsi,t,axis=0)
    Gamman = np.real(Cold - 0.25 * Cold @ np.conj(Mn.T) @ iSigmay @ Mn @ Cold) 
    
    iC = np.linalg.inv(C)
    g = np.real( np.linalg.slogdet(C)[1] + np.dot(np.conj(U), iC @ U) + np.trace(iC @ Gamman) )
                    
    return g