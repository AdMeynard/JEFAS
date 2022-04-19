function C = calc_cov_synth_sharp(M_psi,S,omega,scales,theta)

omega0 = pi/2 ; % the wavelet mode is at Fs/4
vtheta = sqrt( interp1(omega,S,2^(-theta)*omega0./scales,'linear',0) ) ;
C = (vtheta'*vtheta) .* M_psi ;