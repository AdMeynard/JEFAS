function [C,dC] = calc_dcov_sharp(M_psi,S,omega,scales,T,theta)

omega0 = pi/2 ; % the wavelet mode is at Fs/4
vtheta = sqrt( interp1(omega,S,2^(theta)*omega0*scales,'linear',0) ) ;
C = (vtheta'*vtheta) .* M_psi ;

dS = [0 diff(S)]*T ;
dStheta = interp1(omega,dS,2^(theta)*omega0*scales,'linear',0) ;
dvtheta = log(2) * (2^(theta)*scales.*dStheta)./(2*vtheta) ;
tmp = dvtheta'*vtheta ;
dC = (tmp+tmp') .* M_psi ;