function [C,dC] = calc_dcov_sharp(M_psi,S,omega,scales,T,theta)
% CALC_DCOV_SHARP Computation the sharp covariance matrix
% usage:	C = calc_cov_synth_sharp(M_psi,S,omega,scales,theta)
%
% Input:
%   M_psi: matrix output of the function BAS_CALC_COV
%   S: spectrum of X (row vector)
%   omega: vector of frequencies
%   scales: vector of scales
%   theta: Time warping parameter 
%
% Output:
%   C: covariance matrix
%   dC: derivative of the covariance matrix w.r.t. theta

omega0 = pi/2 ; % the wavelet mode is at Fs/4
vtheta = sqrt( interp1(omega,S,2^(theta)*omega0*scales,'linear',0) ) ;
C = (vtheta'*vtheta) .* M_psi ;

dS = [0 diff(S)]*T ;
dStheta = interp1(omega,dS,2^(theta)*omega0*scales,'linear',0) ;
dvtheta = log(2) * (2^(theta)*scales.*dStheta)./(2*vtheta) ;
tmp = dvtheta'*vtheta ;
dC = (tmp+tmp') .* M_psi ;