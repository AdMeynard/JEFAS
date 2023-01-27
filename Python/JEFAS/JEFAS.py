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
import cwtJEFAS
import scipy.optimize as scpo

#%%

def JEFAS(y,Dt,wavparam,freqmin,freqmax,Nf,wavtyp = 'sharp',AM = False,activityThres = 0.05,Nit = 10,errtol = 5e-3,disp = 0 ,**kargs):
    """ JEFAS:	Alternate estimation of the deformations and the spectrum (JEFAS)
    usage:	aML,dgammaML, Sx, crit = JEFAS(y,Dt,wavparam,freqmin,freqmax,Nf,wavtyp = 'sharp',AM = False,activityThres = 0.05,NsS = 110,NfS = 2500,Nit = 10,errtol = 5e-3,disp = 0 ,**kargs)
    
    Parameters
    ----------
    y : array_like
        Signal to analyze
    Dt : int
        Subsampling step for estimations (>=1)
    wavparam : float 
        Wavelet parameter for the time warping estimation (cf. cwtJEFAS)
    freqmin : float
        Minimum normalized frequency for the CWT computing
    freqmax :  float
        Maximum normalized frequency for the CWT computing
    Nf : int
        Number of frequencies/scales for the CWT computing
    wavtyp : string, optional
        Wavelet type (cf. cwtJEFAS), default is wavtyp='sharp'
    AM : bool, default: False
        Whether or not the amplitude modulation is estimated
    activityThres : int, optional
        Activity threshold for JEFAS estimation (0<=activityThres<1), default is 0.05
    Nit : int, optional
        Maximum number of iterations of the alternate algorithm (default: Nit = 10)
    errtol : float, optional
        Relative error in TW, AM and spectrum estimates between iterations that is acceptable for convergence (default is errtol = 0.5%)
    disp : int, optional
        If disp = 0, disable relative update evolution (default: disp = 1)
    
    Other Parameters
    ----------
    TWinit : array_like, optional
        Initial estimation of gamma'(t), default is TWinit[t]=1 for all t
    iterationsTW : int, optional 
        Maximum number of iterations of the optimization algorithm (default: itTW = 6)
    stopTW :  float, optional
        Asbsolute error in the TW parameter between iterations of te optimization algorithm (default: stopTW = 2e-2)
    AMinit : array_like, optional
        Initial estimation of a(t), default is AMinit[t]=1 for all t
    wavparamAM : float, optional
        Wavelet parameter for AM and spectrum estimations (cf. cwtJEFAS), default is wavparamAM=10*wavparam
    regularizationAM : 
        Regularization parameter of the covariance matrix in the AM estimation
    NsS : int, optional
        Number of scales for the wavelet-based spectrum estimator, default is NsS=100
    NfS : int, optional 
        Number of frequencies where the spectrum is estimated (default: Nf = 2500)
      
    Returns
    ----------
    aML : ndarray
        Estimation of the amplitude modulation function, only if AM=True
    dgammaML : ndarray
        Estimation of the time warping function
    Sx : ndarray
        Estimation of the spectrum of the underlying stationary signal x
    crit : ndarray
        Evolution of the relative error in TW, AM and spectrum estimates
    """
      
    #%% Time warping parameters and initialization
    
    T = len(y)
    
    scalesTW = (freqmax/freqmin)**(np.arange(Nf)/(Nf-1)) / (4*freqmax)

    WyTW = cwtJEFAS.cwt(y,scalesTW,wavtyp,wavparam) # Wavelet transform for thetaTW estimation
    
    paramTW = {}
    if 'iterationsTW' in kargs:
        paramTW['iterations'] = kargs['iterationsTW']
    else:
        paramTW['iterations'] = 6
    if 'stopTW' in kargs:
        paramTW['stop'] = kargs['stopTW']
    else:
        paramTW['stop'] = 2e-2
        
    if 'TWinit' in kargs:
        dgamma0 = kargs['TWinit']
    else:
        dgamma0 = np.ones(T)
    
    #%% Amplitude modulation parameters
        
    if 'wavparamAM' in kargs:
        wavparamAM = kargs['wavparamAM']
    else:
        wavparamAM = 10*wavparam   
    
    if AM:
        if 'NfAM' in kargs:
            NfAM = kargs['NfAM']
        else:
            NfAM = 3*Nf
        if 'regularizationAM' in  kargs:
            r = kargs['regularizationAM']
        else:
            r = 1e-5
        
        scalesAM = (freqmax/freqmin)**(np.arange(NfAM)/(NfAM-1)) / (4*freqmax)
        WyAM = cwtJEFAS.cwt(y,scalesAM,wavtyp,wavparamAM) # Wavelet transform for thetaAM estimation
        
    if 'AMinit' in kargs:
        a0 = kargs['AMinit']
    else:
        a0 = np.ones(T)
        
    if (len(dgamma0)!=T) or (len(a0)!=T):
        raise NameError('The initial deformation functions must be of the same length as the signal.')
    
    #%% Spectrum parameters
    if 'NfS' in kargs:
        NfS = kargs['NfS']
    else:
        NfS = 2500
    if 'NsS' in  kargs:
        NsS = kargs['NsS']
    else:
        NsS = 100
    
    scalesS = (0.5/freqmin)**(np.arange(NsS)/(NsS-1)) / 2
        
    WyS = cwtJEFAS.cwt(y,scalesS,wavtyp,wavparamAM) # Wavelet transform for Sx estimation
    
    sigmax = np.var(y)
    
    #%% Active time instants and initializations
    
    act = sigDetection(y,activityThres,scalesTW,wavtyp,wavparamAM)
    timesAct, = np.where( np.diff(act)==1 )
    
    times = np.arange(0,T,Dt)
    times = np.append(times[act[times]==1], timesAct)
    times = np.unique(times)
    
    thetaTW = np.log2(dgamma0[times]) # Initialize thetaTW
    thetaAM = a0[times]**2 # Initialize thetaAM
    Sx = estimSpectrum(WyS,scalesS,times,thetaTW,thetaAM,NfS,sigmax) # initialize Sx
    
    #%% Alternate algorithm
        
    n = 1
    T_est = len(thetaAM)
    tm = int(0.05*T_est) # prevent edge effect from acting on convergence
    tM = int(0.95*T_est) # prevent edge effect from acting on convergence
    errTW = np.inf
    errAM = np.inf
    crit = []
    while (n<=Nit) and ((errTW>errtol) or (errAM>errtol)):
        
        # Step 1: Time warping estimation
        thetaTW_old = thetaTW
        thetaTW0 = thetaTW_old[0]
        thetaTW = estimTW(thetaTW0,thetaAM,times,act,WyTW,scalesTW,Sx,wavtyp,wavparam,paramTW,Dt)
        
        # Step 2: AM estimation
        thetaAM_old = thetaAM
        if AM:
            thetaAM = estimAM(thetaTW,times,act,WyAM,Sx,r,scalesAM,wavtyp,wavparamAM,Dt)
          
        # Step 3: Spectrum estimation
        Sx_old = Sx
        Sx = estimSpectrum(WyS,scalesS,times,thetaTW,thetaAM,NfS,sigmax)
        
        # Relative error evolution
        if np.any(thetaTW_old):
            errTW = np.sum((thetaTW[tm:tM] - thetaTW_old[tm:tM])**2)/np.sum(thetaTW_old[tm:tM]**2)
        else:
            errTW = np.inf
        errAM = np.sum((thetaAM[tm:tM] - thetaAM_old[tm:tM])**2)/np.sum(thetaAM_old[tm:tM]**2)
        errS = np.sum((Sx - Sx_old)**2)/np.sum(Sx_old**2)
        crit = np.append(crit, [errTW, errAM, errS])
        if disp:
            print(f"Iteration {n}")
            print(f"Relative update TW: {100*errTW:.2f} %")
            if AM:
                print(f"Relative update AM: {100*errAM:.2f} %")
            print("")
    
        n = n+1
    
    aML = np.interp(np.arange(T),times,np.sqrt(np.abs(thetaAM)),left=1,right=1)*act + (1-act) # interpolate where the signal is active
    dgammaML = np.interp(np.arange(T),times,2**thetaTW,left=1,right=1)*act + (1-act)
    
    crit = np.reshape(crit,(n-1,3))
    
    if AM:
        return aML,dgammaML, Sx, crit
    else:
        crit = crit[:,(0,2)]
        return dgammaML, Sx, crit

#%%

def sigDetection(y,activityThres,scales,wavtyp,wavparam):
    """SIGDETECTION Detect the activity of the signal
    usage:	act = sigDetection(y,activityThres,scales,wavtyp,wavparam)
    
    Input:
      y: signal to analyze
      activityThres: threshold do detect activity
      scales: vector of scales for wavelet transform
      wavtyp: cf. cwtJEFAS
      wavparam: cf. cwtJEFAS
    
    Output:
      act: is one when the signal is active, zero elsewhere
      """
    
    T = len(y)
    act = np.ones(T)
    
    Wy = cwtJEFAS.cwt(y,scales,wavtyp,wavparam)
    evolPower = np.sum( np.abs(Wy)**2, axis = 0 )
    
    medPower = np.median(evolPower)
    threshold = activityThres*medPower
    
    act = np.where(evolPower > threshold, act, 0)
    return act

#%%

def estimSpectrum(Wy,scales,times,thetaTW,thetaAM,Nf,sigmax):
    """ESTIMSPECTRUM:  Spectrum estimation
    usage:	hatSx = stat_from_cwt(Wy,scales,Dt,thetaTW,thetaAM,Nf,sigmax)
    
    Input:
      Wy: wavelet transform of the signal
      scales: scales on which the CWT is calculated
      times: vector of time for theta estimation
      thetaTW: time warping estimator
      thetaAM: Amplitude modulation estimator
      Nf : number of discrete frequencies
      sigmax : signal power
    
    Output:
      hatSx : estimated spectrum for (2*Nf-1) frequencies between 0 and Fs
      """

    N = len(thetaTW)
    Ms,T = np.shape(Wy)
    
    scales_log = np.log2(scales)
    
    mWy = np.abs(Wy)**2 # scalogram of y
    
    mWx = np.zeros((Ms,N))
    n = 0
    for t in times:
        scalesX = scales_log + thetaTW[n]
        mWx[:,n] = np.interp(scales_log,scalesX,mWy[:,t],left=0,right=0)
        n = n+1

    mWx = mWx/np.sqrt(thetaAM)
    
    weight = N - np.sum(mWx==0,axis=1)
    hatSx = np.sum(mWx,axis=1)/(weight) # time mean of the scalogram of x
    hatSx = np.flip(np.append(hatSx,0)) # add zero frequency Sx(0)=0
    
    xi0 = 1/4 # normalized central frequency of |\hat\psi|^2
    freqs = xi0/scales
    freqs = np.flip(np.append(freqs, 0)) # add zero frequency
    
    freq = np.linspace(0,1/2,Nf) # Nf freq from 0 to Fs/2
    hatSx = np.interp(freq,freqs,hatSx,left=0,right=0)
    hatSx = np.append( hatSx, np.flip(hatSx[1:]) ) # real signal
    
    hatSx = 2*sigmax*Nf*hatSx/np.sum(hatSx) # normalization
    hatSx = np.where(np.isnan(hatSx),0, hatSx)
    
    return hatSx

#%%

def estimTW(thetaTW0,thetaAM,times,act,Wy,scales,S,wavtyp,wavparam,paramTW,Dt):
    """ESTIMTW	Maximum likelihood estimation of the time warping parameter
    usage:	thetaTW = estimTW(thetaAM,Wy,scales,S,wavtyp,wavparam,paramTW,Dt)
    
    Input:
      thetaTW0: initial guess of the TW parameter at time t=0
      thetaAM: guess of AM parameter
      times: vector of times for theta estimation
      act: instants where the signal is active
      Wy: wavelet transform of the signal
      scales: scales on which the CWT is calculated
      S: current guess of the spectrum of z
      wavtyp: type of wavelet (cf. cwt)
      wavparam: cf. cwt
      paramTW: dictionnary of 3 terms, namely
          itTW: maximum number of iteration for each gradient descent
          stopTW: stopping criterion for the gradient descent
      Dt: subsampling rate used for the estimation of thetaTW
    
    Output:
      thetaTW : estimation of the time warping operator
      """

    itTW = paramTW['iterations']
    stopTW = paramTW['stop']

    _,Tf = np.shape(Wy)
    T_hor = len(times)
    thetaTW = np.zeros(T_hor)
    
    Nf = len(S)
    Mpsi = bas_calcCov(scales,wavtyp,wavparam,Nf) # matrices for covariance computations 
    
    theta = thetaTW0
    
    
    # At each time, we run an optimization method
    for n in range(T_hor):
        t = times[n]
        S_AM = thetaAM[n]*S
        
        inst = np.arange(np.maximum(0,(t+1-np.ceil(Dt/2))), np.minimum((t+np.floor(Dt/2))+1,Tf),dtype=int)
        inst = inst[act[inst]==1]
        U = Wy[:,inst] # column(s) we are interested in (we compute a mean over Dt columns)
        
        theta = scpo.fmin(loglh,theta,args=(U,Mpsi,S_AM),xtol=stopTW,maxiter=itTW,disp=0)
        thetaTW[n] = theta
        
    thetaTW = np.log2(2**thetaTW/np.mean(2**thetaTW)) # <gamma'>=1
    
    return thetaTW

#%%

def estimAM(thetaTW,times,act,Wy,S,r,scales,wavtyp,wavparam,Dt):
    """ESTIMAM	Maximum likelihood estimation of the AM parameter
    usage:	thetaAM = estimAM(thetaTW,Wy,S,r,scales,wavtyp,wavparam,Dt)
    
    Input:
      thetaTW: guess of the time warping parameter
      times: vector of times for theta estimation
      act: instants where the signal is active
      Wy: wavelet transform of the signal
      S: current guess of the spectrum
      r: regularization parameter
      scales: scales on which the CWT is calculated
      wavtyp: type of wavelet (cf. cwt)
      wavparam: cf. cwt
      Dt: subsampling rate used for the estimation of theta_MV
    
    Output:
      thetaAM : estimation of the AM parameter
      """

    Ms,Tf = np.shape(Wy)
    N = len(thetaTW)
    
    Nf = len(S)
    Mpsi = bas_calcCov(scales,wavtyp,wavparam,Nf)
    n = 0
    thetaAM = np.zeros(N)
    for t in times:
        inst = np.arange( np.maximum(0,(t+1-Dt/2)), np.minimum((t+Dt/2)+1,Tf),dtype=int )
        inst = inst[act[inst]==1]
        U = Wy[:,inst]
        _,Mu = np.shape(U)
        C = calcCov(Mpsi,S,thetaTW[n])
        C = (1-r)*C + r*np.eye(Ms) # regularization
        thetaAM[n] = np.real( np.trace(np.conj(np.transpose(U)) @ (np.linalg.inv(C) @ U)) ) / Mu
        n = n+1
    
    thetaAM = thetaAM/Ms
    thetaAM = thetaAM/np.mean(thetaAM[10:-10]) # <thetaAM>=1 avoiding edge effect
    
    return thetaAM

#%%

def loglh(theta,U,Mpsi,S):
    """DLOGLH	Computation of the value of the opposite of the loglikelihood
    usage:	g = loglh(theta,U,Mpsi,M_tmpdpsi,S)
    
    Input:
      theta: value of the parameter
      U: data vector (column of the wavelet transform)
      Mpsi: first matrix output of the function BAS_CALC_DCOV
      S: current guess of the spectrum of z
    
    Output:
      g : value of the  the opposite of the loglikelihood
      """


    _,M = np.shape(U)
    
    C = calcCov(Mpsi,S,theta) # estimated covariance matrix and its derivative
    v = np.linalg.inv(C) @ U
    g = np.linalg.slogdet(C)[1] + np.trace(np.conj(np.transpose(U)) @ v).real / M
                                         
    return g

#%%

def calcCov(Mpsi,S,theta):

    T = len(S)
    omega = np.arange(T)*2*np.pi/T
    
    Stheta = np.interp(2**(-theta)*omega,omega,S,left=S[-1],right=S[-1])
    U = np.sqrt(Stheta) * Mpsi
    C = (U @ np.conj(np.transpose(U))) / T
         
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

# #%%

# def calc_dcov(Mpsi,M_tmpdpsi,S,theta):

#     T = len(S)
#     omega = np.arange(T)*2*np.pi/T
    
#     Stheta = np.interp(2**(-theta)*omega,omega,S,left=S[-1],right=S[-1])
#     U = np.sqrt(Stheta) * Mpsi
#     V = np.sqrt(Stheta) * M_tmpdpsi
    
#     UUt = U @ np.conj(np.transpose(U))
#     UVt = U @ np.conj(np.transpose(V))
    
#     C = (UUt)/T
#     dC = np.log(2)*(UUt + UVt + np.conj(np.transpose(UVt)))/T
    
#     return C,dC


# def bas_calc_dcov(scales,wavtyp,wavparam,T):

#     omega = np.arange(T)*2*np.pi/T
    
#     M_tmp = np.outer(scales,omega)
#     scales = np.reshape(scales,(len(scales),1))
    
#     if wavtyp=='sharp':
#         eps = np.finfo(np.float64).eps # numpy precision
#         par = -2*wavparam;
#         M_ftmp1 = (np.pi/(2*M_tmp+eps))
#         M_ftmp2 = (2*M_tmp/np.pi)
#         Mpsi = np.exp( par*(M_ftmp1 + M_ftmp2 - 2) )
#         M_tmpdpsi = par*( M_ftmp2 - M_ftmp1 ) * Mpsi
            
#     elif  wavtyp=='dgauss':
#         Cst = 4*wavparam/(np.pi**2)
#         K = (2/np.pi)**wavparam*np.exp(wavparam/2)
#         M_tmp_pow = M_tmp**wavparam
#         Mpsi = K * M_tmp_pow * np.exp(-Cst*M_tmp**2/2)
#         M_tmpdpsi = (wavparam - Cst*M_tmp**2) * Mpsi
    
#     Mpsi = np.sqrt(scales) * Mpsi
#     M_tmpdpsi = np.sqrt(scales) * M_tmpdpsi

#     return Mpsi,M_tmpdpsi

