function y = warp(x,gamma,dgamma)
%INV_WARP  time warp the signal by deformation function gamma
% usage:	y = warp(x,gamma,dgamma)
%
% Input:
%   y: signal to unwarp
%   gamma: time warping function
%   dgamma: derivative of gamma
% 
% Output:
%   x : x time warped by gamma

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
% Created: 2018-24-05

N = length(x);
ts = linspace(0,1,N);

y = sqrt(dgamma).*interp1(ts,x,gamma,'spline');
end