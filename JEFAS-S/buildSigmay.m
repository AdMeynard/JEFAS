function Sigmay = buildSigmay(MMSigmay,T,TT,Delta)
%BASELINEWARPEST Concatenate short covariance matrices of subsignals to
%create the full covariance matrix of the signal
% usage:	Sigmay = buildSigmay(MMSigmay,T,TT,Delta)
%
% Input:
%   MMSigmay :
%   T :
%   TT :
%   Delta :
% 
% Output:
%   Sigmay : 

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
% Created: 2019-10-01


K = length(MMSigmay);

for k = 1:K % k-th time block
    nn = ( (k-1)*Delta+1 ):( (k-1)*Delta+TT );
    nn = nn((nn>=1)&(nn<=T)); % instants of the k-th block
    NN = length(nn);
    
    Sigmay(nn,nn) = MMSigmay{k} ;
end
