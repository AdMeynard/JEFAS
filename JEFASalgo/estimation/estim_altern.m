function [aML,dgammaML, Sx, crit] = estim_altern(y,Dt,ratio,dgamma0,a0,paramWAV,paramWP,paramAM,paramS,stop_crit,varargin)
%ESTIM_ALTERN	Alternate estimation of the deformations and the spectrum (JEFAS)
% usage:	[aML,dgammaML, Sx, crit] = estim_altern(y,Dt,dgamma0,a0,paramWAV,paramWP,paramAM,paramS,stop_crit,Nit,disp)
%
% Input:
%   y: signal to analyze
%   Dt: subsampling step for estimations (integer >=1)
%   ratio: activity threshold for JEFAS estimation
%   dgamma0: initial estimation of gamma'(t)
%   a0: initial estimation of a(t)
%   paramWAV: cell of three entries: {wav_typ,wav_param,wav_paramWP} where
%       wav_typ: wavelet type (cf. cwt_JEFAS)
%       wav_param: wavelet parameter for AM and spectrum estimations (cf. cwt_JEFAS)
%       wav_paramWP: wavelet parameter for time warping estimation(cf. cwt_JEFAS)
%   paramWP: cell of three entries: {scalesWP,itWP,stopWP} where
%       scalesWP: vector of scales for the WP estimation
%       dmaxWP (optional): maximum accepted difference between two consecutive terms of thetaWP
%       itWP (optional): maximum number of gradient iterations per instant (default: itWP = 6)
%       stopWP (optional): minimum gradient innovation (default: stopWP = 2e-2)
%   paramAM: cell of 1 to 3 entries: paramAM = {AMopt,scalesAM,r} where
%       AMopt: if AM is not estimated AMopt='no AM' => paramAM = {'no AM'}. Esle AMopt='AM' and:
%       scalesAM: vector of scales for the AM estimation
%       r: regularization parameter of the covariance matrix in the AM estimation
%   paramS: cell of 2 entries: paramS = {scalesS,Nf} where
%       scalesS: vector of scales for the spectrum estimation
%       Nf (optional): number of frequencies where the spectrum is estimated (default: Nf = 2500)
%   stop_crit: stopping criterion for the alternate estimation
%   Nit (optional): maximum number of iterations of the alternate algorithm (default: Nit = 10)
%   disp (optional): if disp = 0 disable relative update evolution (default: disp = 1)
% 
% Output:
%   aML : estimation of the amplitude modulation function
%   dgammaML : estimation of the time warping function
%   Sx : Estimation of the spectrum of the underlying stationary signal x
%   crit : evolution of the stopping criterion

% Copyright (C) 2017 Adrien MEYNARD
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Author: Adrien MEYNARD
% Created: 2017-12-19

%% Some parameters
T = length(y);
if (length(dgamma0)~=T)||(length(a0)~=T)
    error('The initial deformation functions must have the same length as the signal.')
end
wav_typ = paramWAV{1};
wav_param = paramWAV{2};
wav_paramWP = paramWAV{3};

%% Time warping parameters
scalesWP = paramWP{1};
WyWP = cwt_JEFAS(y,scalesWP,wav_typ,wav_paramWP); % Wavelet transform for thetaWP estimation

if length(paramWP)==1
    dmaxWP = 0.5;
    itWP = 6;
    stopWP = 2e-2;
elseif length(paramWP)==2
    dmaxWP = paramWP{2};
    itWP = 6;
    stopWP = 2e-2;
elseif length(paramWP)==3
    dmaxWP = paramWP{2};
    itWP = paramWP{3};
    stopWP = 2e-2;
elseif length(paramWP)==4
    dmaxWP = paramWP{2};
    itWP = paramWP{3};
    stopWP = paramWP{4};
else
    error('paramWP must contain 1, 2 or 3 entries')
end

%% Amplitude modulation parameters
AMopt = paramAM{1};

if strcmpi(AMopt,'AM')
    scalesAM = paramAM{2};
    WyAM = cwt_JEFAS(y,scalesAM,wav_typ,wav_param); % Wavelet transform for thetaAM estimation
    r = paramAM{3};
    
elseif ~strcmpi(AMopt,'no AM')
    error('The first entry of paramAM must be ''no AM'' or ''AM''.');
end

%% Spectrum parameters
scalesS = paramS{1};
WyS = cwt_JEFAS(y,scalesS,wav_typ,wav_param); % Wavelet transform for Sx estimation

if length(paramS)==2
    Nf = paramS{2};
else
    Nf = 2500;
end

sigmax = var(y);

%% Active time instants and initializations
act = sigDetection(y,ratio,scalesWP,wav_typ,wav_param);
timesAct = find( subplus( diff(act) ) );

times = 1:Dt:T;
times = [times(act(times)==1) timesAct];
times = sort(times);

thetaWP = log2(dgamma0(times)); % Initialize thetaWP
thetaAM = a0(times).^2; % Initialize thetaAM
Sx = estim_spectrum(WyS,scalesS,times,thetaWP,thetaAM,Nf,sigmax); % initialize Sx

%% Alternate algorithm
if isempty(varargin)
    Nit = 10;
else
    Nit = varargin{1};
end

if (nargin==12)&&(varargin{2}==0)
    disp = 0; % disable evolution criterion display
else
    disp = 1; % enable evolution criterion display
end

n = 1;
T_est = length(thetaAM);
tm = ceil(0.05*T_est); % prevent edge effect from acting on convergence
tM = floor(0.95*T_est); % prevent edge effect from acting on convergence
errWP = inf;
errAM = inf;
crit = [];
while (n<=Nit)&&((errWP>stop_crit)||(errAM>stop_crit))
    
    % Step 1: Time warping estimation
    thetaWP_old = thetaWP;
    thetaWP0 = thetaWP_old(1);
    thetaWP = estim_WP(thetaWP0,thetaAM,times,act,WyWP,scalesWP,Sx,wav_typ,wav_paramWP,itWP,stopWP,Dt,dmaxWP);
    
    % Step 2: AM estimation
    thetaAM_old = thetaAM;
    if strcmpi(AMopt,'AM')
        thetaAM = estim_AM(thetaWP,times,act,WyAM,Sx,r,scalesAM,wav_typ,wav_param,Dt);
    end
      
    % Step 3: Spectrum estimation
    Sx_old = Sx;
    Sx = estim_spectrum(WyS,scalesS,times,thetaWP,thetaAM,Nf,sigmax);
    
    % Stopping criterion evolution
    errWP = sum((thetaWP(tm:tM) - thetaWP_old(tm:tM)).^2)/sum(thetaWP_old(tm:tM).^2);
    errAM = sum((thetaAM(tm:tM) - thetaAM_old(tm:tM)).^2)/sum(thetaAM_old(tm:tM).^2);
    errS = sum((Sx - Sx_old).^2)/sum(Sx_old.^2);
    crit = [crit; [errWP errAM errS]];
    if disp
        if strcmpi(AMopt,'AM')
            fprintf(' Iteration %i \n Relative update WP: %.2f %% \n Relative update AM: %.2f %%\n\n',[n 100*errWP 100*errAM]);
        elseif strcmpi(AMopt,'no AM')
            fprintf(' Iteration %i \n Relative update WP: %.2f %% \n\n',[n 100*errWP]);
        end
    end

    n = n+1;
end

aML = interp1(times,sqrt(abs(thetaAM)),1:T,'linear',1);
dgammaML = interp1(times,2.^thetaWP,1:T,'linear',1);
end