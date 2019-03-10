function [heapBoptim,dgamma,a,Sx,SIRupdate] = altern_bss_deform_nonstat(z,heapB0,vectau0,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,stopSIR,varargin)
%ALTERN_BSS_DEFORM_NONSTAT	Alternate BSS and estimation of the deformations and the spectrum of the sources
% usage:	[heapBoptim,dgamma,Sx,a,SIRupdate] = altern_bss_deform_nonstat(z,heapB0,vectau0,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,stopSIR,Nit,itMAX,stop0)
%
% Input:
%   z: Observed signals (mixture)
%   heapB0: Initial unmixing matrices (thrid dimension for time)
%   vectau0: vector of times where the initial unmixing matrices of given
%   dgamma0: initial estimations of gamma'(t)
%   Dt: subsampling step for time warping estimations (integer >=1)
%   paramBSS: cell of four entries: {scales,vectau,eps_bss,rBSS}
%       scales: vector of scales for the BSS estimation
%       vectau: vector of times where the mixing matrices is estimated in JEFAS-BSS
%       L_bss: number of iterations in the Newton descent
%   paramWP: cell of three entries: {scalesWP,itWP,stopWP} where
%       scalesWP: vector of scales for the WP estimation
%       itWP (optional): maximum number of gradient iterations per instant (default: itWP = 6)
%       stopWP (optional): minimum gradient innovation (default: stopWP = 2e-2)
%   paramAM: cell of 1 to 3 entries: paramAM = {AMopt,scalesAM,r} where
%       AMopt: if AM is not estimated AMopt='no AM' => paramAM = {'no AM'}. Esle AMopt='AM' and:
%       scalesAM: vector of scales for the AM estimation
%       r: regularization parameter of the covariance matrix in the AM estimation
%   paramS: cell of 2 entries: paramS = {scalesS,Nf} where
%       scalesS: vector of scales for the spectrum estimation
%       Nf (optional): number of frequencies where the spectrum is estimated (default: Nf = 2500)
%   stop_crit: stopping criterion for JEFAS
%   stop_SIR: stopping criterion for JEFAS-BSS
%   Nit (optionnal): maximum number of iterations of the alternate algorithm JEFAS-BSS (default: Nit = 15)
%   itMAX (optional): maximum number of iterations in fmincon for unmixing matrix estimation (default: itMAX = 100)
%   stop0 (optional): stopping tolerance in fmincon for unmixing matrix estimation (default: stop0 = 1e-5)
% 
% Output:
%   heapBoptim: estimated unmixing matrices
%   dgamma: estimated time warping functions
%   a: estimates amplitude modulation functions
%   Sx: estimated spectra
%   SIRupdate: SIR update evolution (stopping criterion)
%
% Copyright (C) 2018 Adrien MEYNARD
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
% Created: 2018-10-18


%% Parameters

[N,T] = size(z);

wav_typ = paramWAV{1};
wav_param = paramWAV{2};

scalesBSS = paramBSS{1};
vectau = paramBSS{2};
L = paramBSS{3};

ratio = 0; % signals active everytime

% fmincon parameters
% nonlcon = @norm_row;

if nargin == 13
    Nit = varargin{1};
else
    Nit = 10;
end

% Compute the CWT of the observations:
for n=1:N
    Wz(:,:,n) = cwt_JEFAS(z(n,:),scalesBSS,wav_typ,wav_param);
end

if length(paramS)==2
    Nf = paramS{2};
else
    Nf = 2500;
end
M_psi = bas_calc_cov(scalesBSS,wav_typ,wav_param,2*Nf-1);

%% Initializations

heapBoptim = heapB0;
dgamma = dgamma0;
a = ones(N,T);
haty = nonstatunmixing(z,heapBoptim,vectau0);

%% JEFAS-BSS
cv = 1;
SIRit = 0;
figure;
while ((cv<=Nit)&&(mean(SIRit)<=stopSIR))
    
    % WP estimation
    dgammaprec = dgamma;
    aprec = a;
    for n=1:N
        y = haty(n,:); % estimated source n
        [aML,dgammaML, Sxn] = estim_altern(y,Dt,ratio,dgammaprec(n,:),aprec(n,:),paramWAV,paramWP,paramAM,paramS,stop_crit,10,0);
        dgamma(n,:) = dgammaML;
        a(n,:) = aML;
        Sx(n,:) = Sxn;
    end
    
    % BSS estimation
    B0 = heapBoptim(:,:,1);
%     heapBoptim = estim_mixingmatrix_nonstat(B0,vectau,Wz,Sx,dgamma,M_psi,epsBSS,nonlcon,options);
    A0 = inv(B0);
    heapBoptim = NEWTON_estim_mixingmatrix_nonstat(A0,Wz,vectau,M_psi,Sx,dgamma,L);
    
    % reordering sources to keep correspondency between dgamma and haty
    oldhaty = haty;
    haty = nonstatunmixing(z,heapBoptim,vectau);
    haty = reordersig(oldhaty,haty);
    
    for n=1:N
        subplot(1,N,n); plot(dgamma(n,:));
    end
    drawnow;
    
    % convergence
    [~,SIRit,~,~] = bss_eval_sources(oldhaty,haty);
    SIRupdate(cv) = sum(SIRit)/N;
    fprintf(' Iteration %i \n SIR update: %.2f dB \n\n',[cv SIRupdate(cv)]);
    
    cv = cv+1;
end
