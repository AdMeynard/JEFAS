function [M_psi,M_tmpdpsi] = bas_calc_dcov_sharp(scales,sigma_sharp)
% BAS_CALC_DCOV_SHARP Computation of the part of the covariance and its derivative
% which are independent of theta
% usage:	[M_psi,M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,T)
%
% Input:
%   scales: vector of scales where the CWT is calculated
% 
% Output:
%   M_psi: matrix necessary for covariance matrix computation
%   M_tmpdpsi: matrix necessary for covariance matrix computation

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


s = log2(scales(:)) ;

M_psi = exp(-(s-s').^2/(2*sigma_sharp^2)) ;
M_tmpdpsi = M_psi ;

end