function thetaAM = estim_AM(thetaWP,times,act,Wy,S,r,scales,wav_typ,wav_param,Dt)
%ESTIM_AM	Maximum likelihood estimation of the AM parameter
% usage:	thetaAM = estim_AM(thetaWP,Wy,S,r,scales,wav_typ,wav_param,Dt)
%
% Input:
%   thetaWP: guess of the time warping parameter
%   times: vector of times for theta estimation
%   act: instants where the signal is active
%   Wy: wavelet transform of the signal
%   S: current guess of the spectrum
%   r: regularization parameter
%   scales: scales on which the CWT is calculated
%   wav_typ: type of wavelet (cf. cwt)
%   wav_param: cf. cwt
%   Dt: subsampling rate used for the estimation of theta_MV
% 
% Output:
%   thetaAM : estimation of the AM parameter

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

[Ms,Tf] = size(Wy);
N = length(thetaWP);

Nf = length(S);
M_psi = bas_calc_cov(scales,wav_typ,wav_param,Nf);
n = 1;
thetaAM = zeros(1,N);
for t=times
    inst = max(1,(t+1-Dt/2)):min((t+Dt/2),Tf);
    inst = inst(act(inst)==1);
    U = Wy(:,inst);
    [~,Mu] = size(U);
    C = calc_cov(M_psi,S,thetaWP(n));
    C = (1-r)*C + r*eye(Ms); % regularization
    thetaAM(n) = real(trace(U'*(C\U)))/Mu;
    n = n+1;
end

thetaAM = thetaAM/Ms;
thetaAM = thetaAM/mean(thetaAM(10:(end-10))); % <thetaAM>=1 avoiding edge effect

end