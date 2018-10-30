function haty = nonstatunmixing(z,heapB,vectau)
% NONSTATUNMIXING Estimate the source order from the estimated unmixing matrices
% usage:	haty = nonstatunmixing(z,heapB,vectau)
%
% Input:
%   z: Observations (mixture)
%   heapB: estimated unmixing matrices (third dimension for time)
%   vectau: vector of times where B is estimated in heapB
%
% Output:
%   haty: estimated sources

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
% Created: 2018-10-30

Sizetau = length(vectau);
T = size(z,2);

% piecewise unmixing
for k = 1:(Sizetau-1)
    tau0 = vectau(k);
    tau1 = vectau(k+1) - 1;
    haty(:,tau0:tau1) = heapB(:,:,k)*z(:,tau0:tau1); % k-th segment
end

tau0 = vectau(Sizetau);
haty(:,tau0:T) = heapB(:,:,Sizetau)*z(:,tau0:T); % last segment

