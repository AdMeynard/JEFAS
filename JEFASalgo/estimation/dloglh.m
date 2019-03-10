function [g,dg] = dloglh(theta,U,M_psi,M_tmpdpsi,S)
% DLOGLH	Computation of the value of the opposite of the loglikelihood and its
% derivative
% usage:	[g,dg] = dloglh(theta,U,M_psi,M_tmpdpsi,S)
%
% Input:
%   theta: value of the parameter
%   U: data vector (column of the wavelet transform)
%   M_psi: first matrix output of the function BAS_CALC_DCOV
%   M_tmpdpsi: second matrix output of the function BAS_CALC_DCOV
%   S: current guess of the spectrum of z
% 
% Output:
%   g : value of the  the opposite of the loglikelihood
%   dg : value of its derivative wrt theta

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

[~,M] = size(U);

switch nargout
    case 1
        C = calc_cov(M_psi,S,theta); % estimated covariance matrix
        g = real(logdet(C)) + real(trace(U'*(C\U)))/M;
        
    case 2
        [C,dC] = calc_dcov(M_psi,M_tmpdpsi,S,theta); % estimated covariance matrix and its derivative
        v = C\U;
        CdC = C\dC;
        g = real(logdet(C)) + real(trace(U'*v))/M;
        dg = trace(CdC) - real(trace(U'*CdC*v))/M;
end

end