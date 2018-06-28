function thetaWP = estim_WP(thetaAM,Wy,scales,S,wav_typ,wav_param,itWP,stopWP,Dt)
%ESTIM_WP	Maximum likelihood estimation of the time warping parameter
% usage:	thetaWP = estim_WP(thetaAM,Wy,scales,S,wav_typ,wav_param,itWP,stopWP,Dt)
%
% Input:
%   thetaAM: guess of AM parameter
%   Wy: wavelet transform of the signal
%   scales: scales on which the CWT is calculated
%   S: current guess of the spectrum of z
%   wav_typ: type of wavelet (cf. cwt)
%   wav_param: cf. cwt
%   itWP: maximum number of iteration for each gradient descent
%   stopWP: stopping criterion for the gradient descent
%   Dt: subsampling rate used for the estimation of thetaWP
% 
% Output:
%   thetaWP : estimation of the warping operator

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

[~,Tf] = size(Wy);
t_hor = 1:Dt:Tf; % estimation instants
T_hor = length(t_hor);
thetaWP = zeros(1,T_hor);

Nf = length(S);
[M_psi, M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,Nf); % matrices for covariance computations 

theta = 0;

options = optimoptions('fmincon','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'MaxIterations',itWP,'StepTolerance',stopWP,'Display','off');
% At each time, we run an optimization method
for n = 1:T_hor
    t = t_hor(n);
    S_AM = thetaAM(n)*S;
    U = Wy(:,max(1,(t+1-Dt/2)):min((t+Dt/2),Tf)); % column(s) we are interested in (we compute a mean over Dt columns)
    
    llh = @(x)dloglh(x,U,M_psi,M_tmpdpsi,S_AM);
    theta = fmincon(llh,theta,[],[],[],[],-1.5,1.5,[],options);
    thetaWP(n) = theta;
end

thetaWP = log2(2.^thetaWP/mean(2.^thetaWP)); % <gamma'>=1
end
