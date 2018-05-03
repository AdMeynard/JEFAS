function z = statAMWP(y,a,dgamma)
%STATAMWP:  estimation of the underlying stationary signal
% usage:	z = stat(y,a,dgamma)
%
% Input:
%   y: non-stationary signal
%   a: amplitude modulation function
%   dgamma: derivative of the time warping function gamma
% Output:
%   z : estimated stationary signal

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

gamma = cumsum(dgamma);
gamma = (gamma - gamma(1))/(gamma(end) - gamma(1)); % normalization such as gamma(end) = 1

Dz = y(:)./(a(:));
z = invwarp(Dz,gamma,dgamma)';

end