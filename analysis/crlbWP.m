function b = crlbWP(thetaWP,thetaAM,S,scalesWP,wav_typ,wav_param)
% CRLBWP: Computation of an estimation of the CRLB
% usage:	b = crlbWP(thetaWP,thetaAM,S,scales,wav_typ,wav_param)
%
% Input:
%   thetaWP: actual value of the time warping parameter
%   thetaAM: actual value of the amplitude modulation parameter
%   S: actual value of the spectrum of X
%   scalesWP: vector of scales for thetaWp estimation
%   wav_typ: type of wavelet
%   wav_param: cf. cwt
%
% Output:
%   b: CRLB wrt thetaWP

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
T = length(S);
b = zeros(1,N);

[M_psi,M_tmpdpsi] = bas_calc_dcov(scalesWP,wav_typ,wav_param,T);
for n = 1:N
    SzAM = thetaAM(n)*S;
    [C,dC] = calc_dcov(M_psi,M_tmpdpsi,SzAM,thetaWP(n));
    b(n) = 2/trace((C\dC)^2); % crlb
end

end