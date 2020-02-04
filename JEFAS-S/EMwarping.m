function [dgammaEST,SxEST, W, nllV] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_param,Dt,TT,Delta,alpha,Nit,thres,itD,stopD,varargin)

T = length(y);
Delta = round(Delta);

[M_psi, M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,T); % To compute the covariance matrices
MatPsi = ifft(M_psi.',[],1);

thetaold = thetaINIT;
if length(varargin) == 1
    theta = varargin{1};
    errINIT = sum( abs( thetaINIT(:) - theta(:) ).^2 );
    fprintf(' Initialization \n Quadratic error: %.3f \n\n', errINIT)
end

options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxIterations',itD,'StepTolerance',stopD,'Display','off'); % fmincon

nbit = 1;
nll = Inf;
stopcrit = Inf;

% EM alternate estimation
while ( (nbit<=Nit) && (stopcrit>thres) )
    
    % Step E: Adapted transform estimation
    x = statAMWP(y,ones(T,1),2.^thetaold);
    Sxest = estim_spec(x,T,alpha); % spectrum
    [W, MMSigmay] = transform_adap(y,sigmay,TT,Delta,M_psi,Sxest,thetaold,MatPsi);
    Sigmay = buildSigmay(MMSigmay,T,TT,Delta); % full signal covariance matrix
    iSigmay = inv(Sigmay); % INVERSION !! DANGER !!
    
    % Stem M: Time-warping estimation
    k = 1;
    for t = 1:Dt:T
        U = W(:,t);
        theta0 = thetaold(t);

        Qemx = @(x)Qem(x, theta0, t, U, M_psi, M_tmpdpsi, Sxest, MatPsi, iSigmay); % Minimized function
        thetaEM(k) = fmincon(Qemx,theta0,[],[],[],[],-0.8,0.8,[],options);
        k = k + 1;
    end
    thetanew = interp1(1:Dt:T,thetaEM,1:T,'linear',thetaEM(end)); % interpolation on all the samples
    
    nllold = nll;
    nll = negloglikelihoodsig(y,Sigmay); % negative log-likelihood (decreasing)
    nllV(nbit) = nll;
    
    stopcrit = nllold - nll;
    
    if length(varargin) == 1
        errEM = sum( abs( thetanew(:) - theta(:) ).^2 );
        fprintf(' Iteration %i \n Quadratic Error: %.3f \n\n', nbit, errEM)
    else
        fprintf(' Iteration %i \n Negative loglikelihood: %.3f \n\n', nbit, nllV(nbit))
    end
    
    nbit = nbit + 1;
    thetaold = thetanew;
end
close; 

thetaEST = interp1(1:Dt:T,thetaEM,1:T,'linear',thetaEM(end));
dgammaEST = 2.^thetaEST;
x = statAMWP(y,ones(T,1),dgammaEST);
SxEST = estim_spec(x,T,alpha); % spectrum
W = transform_adap(y,sigmay,TT,Delta,M_psi,Sxest,thetaold,MatPsi);
