function heapBoptim = NEWTON_estim_mixingmatrix_nonstat(A0,Wz,vectau,M_psi,Sx,dgamma,L)
% ESTIM_MIXINGMATRIX_NONSTAT mixing matrix estimation using Newton
% algorithm
% usage:	heapBoptim = NEWTON_estim_mixingmatrix_nonstat(A0,Wz,vectau,M_psi,Sx,dgamma,L)
%
% Input:
%   A0: initial mixing matrix
%   Wz: CWT of the observations
%   vectau: vector of where unmixing matrix is estimated
%   M_psi: matrix output of the function BAS_CALC_COV
%   Sx: current guess of the spectra of the stationary sources 
%   dgamma: current guess of the time warping functions
%   L: number of iterations of the descent
%
% Output:
%   heapBoptim: estimated unmixing matrices

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
% Created: 2019-03-07

[Ms,~,N] = size(Wz);
Kmat = length(vectau);

heapBoptim = zeros(N,N,Kmat);
a = A0(:); % vectorized mixing matrix
A = A0;

% errmat = eps_bss*ones(N); % prevent brutal variations of the unmixing matrix

delta = 1:N;
k = 1;
for tau0 = vectau
    wz = reshape(Wz(:,tau0,:),Ms,N).';
    % covariance matrices:
    for n=1:N
        C(:,:,n) = calc_cov(M_psi,Sx(n,:),log2(dgamma(n,tau0)));
        M(:,:,n) = real(wz * (C(:,:,n) \ wz') ) ;
    end
    
    % log-likelihood optimization:
    for l=1:L
        aold = a;
        Aold = A;
        
        J = zeros(N^2, N^2); % Jacobian matrix
        for n=1:N
            Mna = M(:,:,n) \ Aold ;
            ga( ((n-1)*N+1) : n*N ) = Aold(:,n).' * Mna - double(delta==n) ;
            J( ((n-1)*N+1) : n*N, ((n-1)*N+1) : n*N ) = Mna.' ;
            J( ((n-1)*N+1) : n*N, : )  = J( ((n-1)*N+1) : n*N, : ) + kron( eye(N), Mna(:,n).' );
        end  
        a = aold - J\ga(:) ; % Newton step
        A = reshape(a,N,N);
    end
    heapBoptim(:,:,k) = inv(A);
    
    k = k+1;
end