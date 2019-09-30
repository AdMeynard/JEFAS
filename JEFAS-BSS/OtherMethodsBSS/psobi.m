function [heapA, heapB] = psobi(z,Ny,vectau)
% PSOBI Piecewise SOBI BSS algorithm
% usage:	[heapA, heapB] = psobi(z,vectau)
%
% Input:
%   z: observed signals
%   vectau: vector of subintervals beginnings
%
% Output:
%   heapA: estimated mixing matrix
%   heapB: estimated unmixing matrix

% Copyright (C) 2018 Adrien MEYNARD
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
% Created: 2018-10-18

[Nz,T] = size(z);
dtau = vectau(2) -vectau(1);
ltau = length(vectau);

heapA = zeros(Nz,Ny,ltau);
heapB = zeros(Ny,Nz,ltau);
k = 1;
for tau=vectau
   Tau = tau:min((tau+dtau-1),T); 
   ztau = z(:,Tau);
   heapA(:,:,k) = sobi(ztau,Ny);
   B = pinv(heapA(:,:,k));
   heapB(:,:,k) = B./sqrt(sum(abs(B).^2,2)); % rows normalization
   k = k+1;
end