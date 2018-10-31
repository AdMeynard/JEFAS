function heapBoptim = estim_mixingmatrix_nonstat(B0,vectau,Wz,Sx,dgamma,M_psi,eps_bss,r,nonlcon,options)
% ESTIM_MIXINGMATRIX_NONSTAT unmixing matrix estimation 
% usage:	heapBoptim = estim_mixingmatrix_nonstat(B0,vectau,Wz,Sx,dgamma,M_psi,eps_bss,r,nonlcon,options)
%
% Input:
%   B0: initial unmixing matrix
%   vectau: vector of where unmixing matrix is estimated
%   Wz: CWT of the observations
%   Sx: current guess of the spectra of the stationary sources 
%   dgamma: current guess of the time warping functions
%   M_psi: matrix output of the function BAS_CALC_COV
%   eps_bss: constant on the bayesian prior
%   r: regularization parameter on the covariance matrix
%   nonlcon: non linear constaint in fmincon for optimization
%   options: options in fmincom for optimization
%
% Output:
% heapBoptim: estimated unmixing matrices

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
Kmat = length(vectau);

heapBoptim = zeros(N,N,Kmat);
Bold = B0;
errmat = eps_bss*ones(N); % prevent brutal variations of the unmixing matrix 
k = 1;
for tau0 = vectau
    % covariance matrices:
    for n=1:N
        Cn = calc_cov(M_psi,Sx(n,:),log2(dgamma(n,tau0)));
        C(:,:,n) = (1-r)*Cn + r*eye(Ms); % regularization
    end
    
    % log-likelihood optimization:
    costB = @(B)cost_lh(B,tau0,Wz,C);
    Boptim = fmincon(costB,Bold,[],[],[],[],Bold-errmat,Bold+errmat,nonlcon,options);
    heapBoptim(:,:,k) = Boptim;
    
    k = k+1;
    Bold = Boptim;
end
