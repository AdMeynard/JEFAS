function y = deformAMWP(z,a,dgamma)
%deformAMWP:  synthesis of a new nonstationary signal
% usage:	y = deformAMWP(z,a,dgamma)
%
% Input:
%   z: stationary signal
%   a: amplitude modulation function
%   dgamma: derivative of the time warping function gamma
% Output:
%   y : synthetic nonstationary signal

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
% Created: 2018-05-24

Td = length(dgamma);
Tz = length(z);
td = linspace(0,1,Td);
tz = linspace(0,1,Tz);
dgamma = interp1(td,dgamma,tz); % matching sizes
a = interp1(td,a,tz); % matching sizes

gamma = cumsum(dgamma);
gamma = (gamma - gamma(1))/(gamma(end) - gamma(1)); % normalization such as gamma(end) = 1

Dz = warp(z,gamma,dgamma)';
y = a(:).*Dz(:);

end