function [g,dg] = Qem(theta, thetaprec, t, U, priorList, M_psi, M_tmpdpsi, S, MatPsi, iSigmay)
%Qem:  evaluation of the function Q to minimize during EM algorithm 
% usage:	[g,dg] = Qem(theta, thetaprec, t, U, priorList, M_psi, M_tmpdpsi, S, MatPsi, iSigmay)
%
% Input:
%   theta: value of the parameter
%   thetaprec: value of the parameter given by the previous EM iteration
%   t: time
%   U: column of the tim-scale representation at time t
%   priorList: shape of the a priori covariance matrix ('wavelet', 'sharp' or 'sparse')
%   M_psi: first matrix output of the function BAS_CALC_DCOV
%   M_tmpdpsi: second matrix output of the function BAS_CALC_DCOV
%   S: spectrum
%   MatPsi: matrix of "wavelets" for synthesis 
%   iSigmay: inverse of the covariance matrix of the signal for synthesis
% Output:
%   g : value of Q(theta, thetaprec)
%   dg : value of dQ(theta, thetaprec)/dtheta

T = length(S);
omega = (0:(T-1))*2*pi/T;

prior = priorList{1} ;
if ~strcmp(prior,'wavelet')
    scales = priorList{2} ;
end

switch nargout
    case 1
        switch prior
            case 'wavelet'
                Cold = calc_cov_synth(M_psi,S,omega,T,thetaprec); % Covariance old
                C = calc_cov_synth(M_psi,S,omega,T,theta); % covariance new
            case 'sharp'
                Cold = calc_cov_synth_sharp(M_psi,S,omega,scales,thetaprec); % Covariance old
                C = calc_cov_synth_sharp(M_psi,S,omega,scales,theta); % covariance new
        end

        Mn = circshift(MatPsi,t-1,1);
        Gamman = real(Cold - 0.25*Cold * Mn' * iSigmay * Mn * Cold) ;

        g = real( logdet(C) + U'*(C\U) + trace(C\Gamman) );
                
    case 2
        switch prior
            case 'wavelet'
                Cold = calc_cov_synth(M_psi,S,omega,T,thetaprec); % Covariance old
                [C,dC] = calc_dcov(M_psi,M_tmpdpsi,S,theta); % estimated covariance matrix and its derivative
            case 'sharp'
                Cold = calc_cov_synth_sharp(M_psi,S,omega,scales,thetaprec); % Covariance old
                [C,dC] = calc_dcov_sharp(M_psi,S,omega,scales,T,theta); % estimated covariance matrix and its derivative 
        end

        Mn = circshift(MatPsi,t-1,1);
        Gamman = real(Cold - 0.25*Cold * Mn' * iSigmay * Mn * Cold) ;
        
        v = C\U;
        CdC = C\dC;
        invCG = C\Gamman;
        g = real(logdet(C)) + real(U'*v) + real(trace(invCG));
        dg = trace(CdC) - real(U'*CdC*v) - real(trace(CdC*invCG));
end

end