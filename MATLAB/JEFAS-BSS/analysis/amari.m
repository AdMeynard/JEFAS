function indexdB = amari(Best, A)
% AMARI evaluate the Amari index from a BSS estimate
% usage:	indexdB = amari(Best, A)
%
% Input:
%   Best: unmixing matrix estimate
%   A: ground truth mixing matrix
%
% Output:
%   indexdB: Amari index in decibel

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

G = Best*A; % should be close to identity
N = size(G,1);

index = (1/( 2*N*(N-1) ))*( sum( sum(abs(G),2)./max(abs(G),[],2) ) + ...
    sum( sum(abs(G),1)./max(abs(G),[],1) ) - 2*N );

indexdB = 10*log10(index);