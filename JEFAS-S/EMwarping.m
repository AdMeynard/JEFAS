function [dgammaEST,SxEST, W, nllV] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_param,priorList,Dt,TT,Delta,TS,alpha,Nit,thres,itD,stopD,varargin)
%EMWARPING	JEFAS-S joint estimation of the time warping function, spectrum, and adaped transform
% usage:	[dgammaEST,SxEST, W, nllV] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_param,priorList,Dt,TT,Delta,TS,alpha,Nit,thres,itD,stopD,theta)
%
% Input:
%   y: signal
%   sigmay: noise variance
%   thetaINIT: initial guess of the WP parameters
%   scales: scales on which the CWT is calculated
%   wav_typ: type of wavelet (cf. cwt)
%   wav_param: wavelet parameter (cf. cwt)
%   priorList: cell of shape {prior, scales} where
%     prior: shape of the a priori covariance matrix ('wavelet', 'sharp' or 'sparse')
%     scales: vector of scales (only when prior~='wavelet' because scales is included in M_psi)
%   Dt: subsampling rate used for the estimation of thetaWP
%   TT: slicing size for W
%   Delta: overlap
%   alpha: width for Welch estimation
%   Nit: maximum number of alternate estimations (iterations)
%   thres: stopping criterion
%   itD: number of iterations for gradient descent
%   stopD: gradient descent: tolerance
%   theta (optionnal): ground-truth value of the WP parameters
% 
% Output:
%   dgammaEST: estimation of the time warping function
%   SxEST: estimation of the underlying spectrum
%   W: time-scale reprenentation
%   nnlV: evolution of the negative log-likelihood through JEFAS-S iterations

T = length(y);
Delta = round(Delta);

[M_psi, ~] = bas_calc_dcov(scales,wav_typ,wav_param,T);
MatPsi = ifft(M_psi.',[],1); % Synthesis coefficients

prior = priorList{1} ;
switch prior % To compute the prior covariance matrices
    case 'wavelet'
        [M_psi, M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,TS) ;
    case 'sharp'
        sigma_sharp = priorList{2} ;
        [M_psi, M_tmpdpsi] = bas_calc_dcov_sharp(scales,sigma_sharp) ;
        priorList{2} = scales ; % sigma_sharp is no longer needed
end

thetaold = thetaINIT;
if length(varargin) == 1
    theta = varargin{1};
    errINIT = sum( abs( thetaINIT(:) - theta(:) ).^2 );
    fprintf(' Initialization \n Quadratic error: %.3f \n\n', errINIT)
end

options.Display = 0 ;
options.MaxIter = itD ;
options.optTol = stopD ;

nbit = 1;
nll = 1/eps ;
errEM = Inf;
stopcrit = Inf;

% EM alternate estimation
while ( (nbit<=Nit) && (stopcrit>thres) )
    
    % Step E: Adapted transform estimation
    x = statAMWP(y,ones(T,1),2.^thetaold); % no AM estimation here
    Sxest = estim_spec(x,TS,alpha); % spectrum
    [W, MMSigmay] = transform_adap(y,sigmay,priorList,TT,Delta,M_psi,Sxest,thetaold,MatPsi);
    Sigmay = buildSigmay(MMSigmay,T,TT,Delta); % full signal covariance matrix
    iSigmay = inv(Sigmay); % INVERSION !! DANGER !!
    
    % Step M: Time-warping estimation
    k = 1;
    for t = 1:Dt:T
        U = W(:,t);
        theta0 = thetaold(t);

        Qemx = @(x)Qem(x, theta0, t, U, priorList, M_psi, M_tmpdpsi, Sxest, MatPsi, iSigmay); % Minimized function
        thetaEM(k) = minFunc(Qemx,theta0,options);
        k = k + 1;
    end
    thetanew = interp1(1:Dt:T,thetaEM,1:T,'linear',thetaEM(end)); % interpolation on all the samples
    
    nllold = nll;
    nll = negloglikelihoodsig(y,Sigmay); % negative log-likelihood (decreasing)
    nllV(nbit) = nll;
    
    stopcrit = nllold - nll ;
    
    if length(varargin) == 1
        errEMold = errEM ;
        errEM = sum( abs( thetanew(:) - theta(:) ).^2 );
        if errEM>errEMold
            break;
        end
        fprintf(' Iteration %i \n Quadratic Error: %.3f \n\n', nbit, errEM)
    else
        fprintf(' Iteration %i \n Negative loglikelihood: %.3f \n\n', nbit, nllV(nbit))
    end
    
    nbit = nbit + 1;
    thetaold = thetanew;
end
close; 

dgammaEST = 2.^thetaold ;
x = statAMWP(y,ones(T,1),dgammaEST);
SxEST = estim_spec(x,T,alpha); % spectrum
W = transform_adap(y,sigmay,priorList,TT,Delta,M_psi,Sxest,thetaold,MatPsi);
