function [heapBoptim,dgamma,a,Sx,errSIR,errSAR,errSDR,estSIR] = altern_bss_deform_nonstat(yv,z,heapB0,vectau0,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,cvmax,stopSIR,options)
%ALTERN_BSS_DEFORM_NONSTAT	Alternate BSS and estimation of the deformations and the spectrum of the sources
% usage:	[heapBoptim,dgamma,Sx,errSIR,errSAR,errSDR,estSIR] = altern_bss_deform_nonstat(yv,z,heapB0,vectau0,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,cvmax,stopSIR,options)
%
% Input:
% 
% Output:
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

[N,T] = size(z);

wav_typ = cell2mat(paramWAV(1));
wav_param = cell2mat(paramWAV(2));

if length(paramS)==2
    Nf = cell2mat(paramS(2));
else
    Nf = 2500;
end

scalesBSS = paramBSS{1};
vectau = paramBSS{2};
epsBSS = paramBSS{3};
rBSS = paramBSS{4};

% Compute the CWT of the observations:
for n=1:N
    Wz(:,:,n) = cwt_JEFAS(z(n,:),scalesBSS,wav_typ,wav_param);
end
M_psi = bas_calc_cov(scalesBSS,wav_typ,wav_param,2*Nf-1);

% Initializations:
heapBoptim = heapB0;
dgamma = dgamma0;
a = ones(N,T);

nonlcon = @norm_row;

haty = nonstatunmixing(z,heapBoptim,vectau0); % revoir cette fonction

[SDR,SIR,SAR,~] = bss_eval_sources(haty,yv); % unmixing quality of the initialization

cv = 1;
SIRit = 0;
figure;
while ((cv<=cvmax)&&(mean(SIRit)<=stopSIR))
    
    % WP estimation
    dgammaprec = dgamma;
    aprec = a;
    for n=1:N
        y = haty(n,:); % estimated source n
        [aML,dgammaML, Sxn] = estim_altern(y,Dt,dgammaprec(n,:),aprec(n,:),paramWAV,paramWP,paramAM,paramS,stop_crit,2,0); % J'ai mis 2 itérations, ça améliore l'algo...
        dgamma(n,:) = dgammaML;
        a(n,:) = aML;
        Sx(n,:) = Sxn;
    end
    
    % BSS estimation
    B0 = heapBoptim(:,:,1);
    heapBoptim = estim_mixingmatrix_nonstat(B0,vectau,Wz,Sx,dgamma,M_psi,epsBSS,rBSS,nonlcon,options);
    
    % reordering sources to keep correspondency between dgamma and haty
    oldhaty = haty;
    haty = nonstatunmixing(z,heapBoptim,vectau);
    haty = reordersig(oldhaty,haty);
    for n=1:N
        Why(:,:,n) = cwt_JEFAS(haty(n,:),scalesBSS,wav_typ,wav_param);
        subplot(1,N,n); imagesc(abs(Why(:,:,n)));
    end
    drawnow;
    
    % convergence ?
    SDR_old = SDR;
    SIR_old = SIR;
    SAR_old = SAR;
    [SDR,SIR,SAR,~] = bss_eval_sources(haty,yv);
    [~,SIRit,~,~] = bss_eval_sources(oldhaty,haty);
    errSDR(cv) = sum(SDR - SDR_old)/N;
    errSIR(cv) = sum(SIR - SIR_old)/N;
    errSAR(cv) = sum(SAR - SAR_old)/N;
    estSIR(cv) = sum(SIRit)/N;
    fprintf(' Iteration %i \n update SDR: %.2f dB \n update SIR: %.2f dB\n update SAR: %.2f dB\n\n',[cv errSDR(cv) errSIR(cv) errSAR(cv)]);
    
    cv = cv+1;
end