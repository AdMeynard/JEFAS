function C = calc_cov_synth(M_psi,S,omega,T,theta)
% CALC_COV_SYNTH Computation the covariance matrix
% usage:	C = calc_cov_synth(M_psi,S,omega,T,theta)
%
% Input:
%   M_psi: matrix output of the function BAS_CALC_COV
%   S: spectrum of X (row vector)
%   omega: vector of frequencies
%   T: signal length
%   theta: Time warping parameter 
%
% Output:
%   C: covariance matrix

Stheta = interp1(omega,S,2^(-theta)*omega,'linear',0);
U = bsxfun(@times,sqrt(Stheta),M_psi);
C = (U*U')/T;