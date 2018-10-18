function A = BSS_QTF(z,pp, eps3, eps4, nn)
% BSS_TFQ BSS using quadratic time-frequency representations
% usage:	A = BSS_TFQ(z,pp, eps3, eps4, nn)
%
% Input:
%   z: observed signals
%   pp:
%   eps3:
%   eps4:
%   nn:
%
% Output:
%   A: estimated mixing matrix

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

N = size(z,1);
[~, Dmat] = selecpts(z, pp, 0, 1, eps3, eps4);

Kd = size(Dmat,3);
for k= 1:Kd
    Dx = real(Dmat(:,:,k));
    [U,~,~] = svd(Dx);
    colA(:,k) = U(:,1);
end

[idx,C] = kmeans(colA',nn); % clustering
nC = histc(idx(:),1:nn); % frequency of each class
for k=1:N
    A(:,k) = C(nC==max(nC),:).';
    nC(nC==max(nC)) = 0;
end