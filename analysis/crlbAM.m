function b = crlbAM(thetaAM,scalesAM)
% CRLBAM: Computation of an estimation of the CRLB
% usage:	b = crlbAM(thetaWP,thetaAM,S,scales,wav_typ,wav_param,noise_param,reg_param)
%
% Input:
%   thetaAM: actual value of the amplitude modulation parameter
%   scales: vector of scales for thetaAM estimation
%
% Output:
%   b: CRLB wrt thetaAM

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

Ms = length(scalesAM);
b = 2*thetaAM.^2/Ms;

end