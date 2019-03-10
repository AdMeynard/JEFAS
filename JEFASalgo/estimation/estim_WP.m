function thetaWP = estim_WP(thetaWP0,thetaAM,times,act,Wy,scales,S,wav_typ,wav_param,itWP,stopWP,Dt,dmaxWP)
%ESTIM_WP	Maximum likelihood estimation of the time warping parameter
% usage:	thetaWP = estim_WP(thetaAM,Wy,scales,S,wav_typ,wav_param,itWP,stopWP,Dt,dmaxWP)
%
% Input:
%   thetaWP0: initial guess of the WP parameter at time t=0
%   thetaAM: guess of AM parameter
%   times: vector of times for theta estimation
%   act: instants where the signal is active
%   Wy: wavelet transform of the signal
%   scales: scales on which the CWT is calculated
%   S: current guess of the spectrum of z
%   wav_typ: type of wavelet (cf. cwt)
%   wav_param: cf. cwt
%   itWP: maximum number of iteration for each gradient descent
%   stopWP: stopping criterion for the gradient descent
%   Dt: subsampling rate used for the estimation of thetaWP
%   dmaxWP: maximum accepted difference between two consecutive terms of thetaWP
% 
% Output:
%   thetaWP : estimation of the time warping operator

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
T_hor = length(times);
thetaWP = zeros(1,T_hor);

Nf = length(S);
[M_psi, M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,Nf); % matrices for covariance computations 

theta = thetaWP0;

options = optimoptions('fmincon','Algorithm','trust-region-reflective','SpecifyObjectiveGradient',true,'MaxIterations',itWP,'StepTolerance',stopWP,'Display','off');
% At each time, we run an optimization method
for n = 1:T_hor
    t = times(n);
    S_AM = thetaAM(n)*S;
    
    inst = max(1,(t+1-ceil(Dt/2))):min((t+floor(Dt/2)),Tf);
    inst = inst(act(inst)==1);
    U = Wy(:,inst); % column(s) we are interested in (we compute a mean over Dt columns)
    
    llh = @(x)dloglh(x,U,M_psi,M_tmpdpsi,S_AM);
    if (t>0.04*Tf)&&(t<0.96*Tf) % outside edges
        theta = fmincon(llh,theta,[],[],[],[],theta-dmaxWP,theta+dmaxWP,[],options);
    else % edge effects
        theta = fmincon(llh,theta,[],[],[],[],-1.5,1.5,[],options);
    end
    thetaWP(n) = theta;
end

thetaWP = log2(2.^thetaWP/mean(2.^thetaWP)); % <gamma'>=1
end
