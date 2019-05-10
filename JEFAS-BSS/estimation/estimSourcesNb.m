function N = estimSourcesNb(z,Delta,thres)
% ESTIMSOURCESNB estimate the sources number
% usage:	Nt = estimSourcesNb(z,Delta,thres)
%
% Input:
%   z: observations
%   Delta: length of the windows
%   thres: threshold for the singular values
%
% Output:
%   N: estimated number of sources

% Copyright (C) 2019 Adrien MEYNARD
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
% Created: 2019-05-10

T = size(z,2);

k = 0;
tf = 0;
while tf<T
    ti= k*Delta+1;
    tf = (k+1)*Delta ;
    tt = ti:tf;
    tt = tt(tt<=T);
    ztt = z(:,tt);
    
    u = svd(ztt);
    Nt(k+1) = length(u(u>thres));
    
    k = k+1;
end

N = mode(Nt);