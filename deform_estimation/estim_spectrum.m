function hatSx = estim_spectrum(Wy,scales,Dt,thetaWP,thetaAM,Nf,sigmax)
%ESTIM_SPECTRUM:  Spectrum estimation
% usage:	hatSx = stat_from_cwt(Wy,scales,Dt,thetaWP,thetaAM,Nf,sigmax)
%
% Input:
%   Wy: wavelet transform of the signal
%   scales: scales on which the CWT is calculated
%   Dt: subsampling step for theta estimation
%   thetaWP: time warping estimator
%   thetaAM: Amplitude modulation estimator
%   Nf : number of discrete frequencies
%   sigmax : signal power
% 
% Output:
%   hatSx : estimated spectrum for (2*Nf-1) frequencies between 0 and Fs

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

N = length(thetaWP);
[Ms,T] = size(Wy);

scales = scales(:);
scales_log = log2(scales);

mWy = abs(Wy).^2; % scalogram of y

mWx = zeros(Ms,N);
n = 1;
for t = 1:Dt:T
    scalesX = scales_log + thetaWP(n);
    mWx(:,n) = interp1(scalesX,mWy(:,t),scales_log,'linear',0);
    n = n+1;
end
thetaAM = thetaAM(:).';
mWx = bsxfun(@rdivide,mWx,sqrt(thetaAM));

weight = N - sum(mWx==0,2);
hatSx = sum(mWx,2)./(weight); % time mean of the scalogram of x
hatSx = [hatSx; 0]; % add zero frequency Sx(0)=0

xi0 = 1/4; % normalized central frequency of |\hat\psi|^2
freqs = xi0./scales;
freqs = [freqs; 0]; % add zero frequency

freq = linspace(0,1/2,Nf); % Nf freq from 0 to Fs/2
hatSx = interp1(freqs,hatSx,freq,'linear',0);
hatSx = [hatSx fliplr(hatSx(2:end))]; % real signal

hatSx = 2*sigmax*Nf*hatSx/sum(hatSx); % normalization
hatSx(isnan(hatSx)) = 0;
end
