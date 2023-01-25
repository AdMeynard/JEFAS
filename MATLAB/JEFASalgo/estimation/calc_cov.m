function C = calc_cov(M_psi,S,theta)
% CALC_COV Computation the covariance matrix
% usage:	C = calc_cov(M_psi,S,theta)
%
% Input:
%   M_psi: matrix output of the function BAS_CALC_COV
%   S: spectrum of X (row vector)
%   theta: Time warping parameter 
%
% Output:
%   C: covariance matrix

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
C = (U*U')/T;

end