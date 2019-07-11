function [g,dg] = Qem(theta, thetaprec, t, U, M_psi, M_tmpdpsi, S, MatPsi, iSigmay)
%Qem:  evaluation of the function Q to minimize during EM algorithm 
% usage:	[g,dg] = Qem(theta, thetaprec, t, U, M_psi, M_tmpdpsi, S, MatPsi, Sigmay)
%
% Input:
%   theta: value of the parameter
%   thetaprec: value of the parameter given by the previous EM iteration
%   M_psi: first matrix output of the function BAS_CALC_DCOV
%   M_tmpdpsi: second matrix output of the function BAS_CALC_DCOV
%   S: spectrum
%   MatPsi: matrix of "wavelets" for synthesis 
%   Sigmay: covariance matrix of the signal for synthesis
% Output:
%   g : value of Q(theta, thetaprec)
%   dg : value of dQ(theta, thetaprec)/dtheta

T = length(S);
omega = (0:(T-1))*2*pi/T;

switch nargout
    case 1
        Cold = calc_cov(M_psi,S,omega,T,thetaprec); % Covariance old
        Mn = circshift(MatPsi,t-1,1);
        rGamman = real(Cold - 0.25*Cold * Mn' * iSigmay * Mn * Cold) ;

        C = calc_cov(M_psi,S,theta); % covariance new
        g = real( logdet(C) + U'*(C\U) + trace(C\rGamman) );
                
    case 2
        Cold = calc_cov(M_psi,S,omega,T,thetaprec); % Covariance old
        Mn = circshift(MatPsi,t-1,1);
        Gammat = Cold - 0.5*Cold * Mn' * iSigmay * Mn * Cold ;
        
        [C,dC] = calc_dcov(M_psi,M_tmpdpsi,S,theta); % estimated covariance matrix and its derivative
        v = C\U;
        CdC = C\dC;
        invCG = C\Gammat;
        g = real(logdet(C)) + real(U'*v) + real(trace(invCG));
        dg = trace(CdC) - real(U'*CdC*v) - real(trace(CdC*invCG));
end

end