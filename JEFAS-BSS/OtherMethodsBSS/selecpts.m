function [ZDmat, Dmat, ptszd, ptsd] = selecpts(z,pp,eps1,eps2,eps3,eps4)
% SELECPTS Points selection for joint zero-diagonalizations (ZD) and
% diagonalizations (D) in the T-F domain
% usage:	[ZDmat, Dmat, ptszd, ptsd] = selecpts(z,pp,eps1,eps2,eps3,eps4)
%
% Input:
%   z: observed signals
%   pp: Time subsampling for QTF representations
%   eps1: imaginary part threshold on points selections for ZD
%   eps2: real part threshold on points selections for ZD
%   eps3: imaginary part threshold on points selections for D
%   eps4: real part threshold on points selections for D
%
% Output:
%   ZDmat: matrices to be zero-diagonalized
%   Dmat: matrices to be diagonalized
%   ptszd: time-frequency points indices where ZD is performed
%   ptsd: time-frequency points indices where D is performed

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

[N,T] = size(z); % original sig
tt = 1:pp:T; % subsampling in order to have correct size matrices
TT = length(tt); % subsampledsig

% RTFQ construction
Dz = zeros(TT,TT,N,N);
for n = 1:N
    for m = 1:N
        Dznm = tfrpwv([hilbert(z(n,:)).', hilbert(z(m,:)).'],tt,TT);
        Dz(:,:,n,m) = Dznm;
    end
end
Mat = permute(Dz,[3 4 1 2]);
Mat = reshape(Mat,N,N,TT^2);

% Points selection
iDz = sum(sum(abs(imag(Dz)).^2,4),3);
rDz = sum(sum(abs(real(Dz)).^2,4),3);
% ZD
ptszd = find((iDz > eps1)&(iDz > eps2*rDz));
ZDmat = Mat(:,:,ptszd);
% D
ptsd = find((iDz < eps3)&(rDz > eps4));
Dmat = Mat(:,:,ptsd);
