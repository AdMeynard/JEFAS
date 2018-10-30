function g = cost_lh(B,tau,Wz,Cn)
% COST_LH log_likelihood evaluation
% usage:	g = cost_lh(B,tau,Wz,Cn)
%
% Input:
% B: unmixing matrix
% tau: evaluation time
% Wz: wavelet transform of the observations
% Cn: covariances matrices for each source
%
% Output:
% g: log-likelihood value

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

[Ms,~,N] = size(Wz);

wz = reshape(Wz(:,tau,:),Ms,N).';
hat_wy = B*wz;

g = -Ms*real(logdet(B));
for n = 1:N       
	C = Cn(:,:,n);
    ing = real(hat_wy(n,:)*(C\(hat_wy(n,:)')));
    g = g + ing;
end