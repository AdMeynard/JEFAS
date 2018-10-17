function haty = association_mariage(haty,vectau,Deltat)
% ASSOCIATION_MARIAGE Piecewise reordering of the sources 
% usage:	haty = association_mariage(haty,vectau,Deltat)
%
% Input:
% haty: unordered sources
% vectau: vector of connection times
% Deltat: number of sample to consider and compare before and after the connexion time
%
% Output:
% haty: reordered sources

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
% Created: 2018-07-11

N = size(haty,1);
K = length(vectau);

for k = 1:K
    tau1 = vectau(k);
    [Prefm, Prefw] = matrice_preference(haty,tau1,Deltat);
    
    stablematch = galeshapley(N,Prefm,Prefw);
    haty(:,1:(tau1-1)) = haty(stablematch,1:(tau1-1)); % reordering
end