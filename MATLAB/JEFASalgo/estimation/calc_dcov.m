function [C,dC] = calc_dcov(M_psi,M_tmpdpsi,S,theta)
% CALC_DCOV Computation the covariance matrix and its derivative
% usage:	[C,dC] = calc_dcov(M_psi,M_tmpdpsi,S,theta)
%
% Input:
%   M_psi: first matrix output of the function BAS_CALC_DCOV
%   M_tmpdpsi: second matrix output of the function BAS_CALC_DCOV
%   S: spectrum of X (row vector)
%   theta: Time warping parameter 
%
% Output:
%   C: covariance matrix
%   dC: derivative of the covariance matrix

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

T = length(S);
omega = (0:(T-1))*2*pi/T;

Stheta = interp1(omega,S,2^(-theta)*omega,'linear',S(end));
U = bsxfun(@times,sqrt(Stheta),M_psi);
V = bsxfun(@times,sqrt(Stheta),M_tmpdpsi);

UUt = U*U';
UVt = U*V';

C = (UUt)/T;
dC = log(2)*(UUt + UVt + UVt')/T;

end