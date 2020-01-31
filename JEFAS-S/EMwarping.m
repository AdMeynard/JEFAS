function [dgammaEST,SxEST, W, nllV] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_param,Dt,TT,Delta,alpha,Nit,thres,itD,stopD,varargin)

T = length(y);
Delta = round(Delta);

[M_psi, M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,T); % Pour calculer les matrices de covariance
% MatPsi = real(ifft(M_psi.',[],1)); % matrice d'ondelettes
MatPsi = ifft(M_psi.',[],1);
% MatPsiD = MatPsi(1:TT,:);

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
while ( (nbit<=Nit) && (stopcrit>thres) )
    
    % Calcul de W: etape E de EM
    x = statAMWP(y,ones(T,1),2.^thetaold);
    Sxest = estim_spec(x,T,alpha); % spectre
    [W, MMSigmay] = transform_adap(y,sigmay,TT,Delta,M_psi,Sxest,thetaold,MatPsi);
    Sigmay = buildSigmay(MMSigmay,T,TT,Delta); % matrice de tout le signal
    iSigmay = inv(Sigmay); % son inverse !! Calcul compliqué !!
    
    % Estimation de thetaEM nouveau: etape M de EM
    k = 1;
    for t = 1:Dt:T
        U = W(:,t); % Wvraie(:,t);
        theta0 = thetaold(t);

        Qemx = @(x)Qem(x, theta0, t, U, M_psi, M_tmpdpsi, Sxest, MatPsi, iSigmay); % fonction a minimiser
        thetaEM(k) = fmincon(Qemx,theta0,[],[],[],[],-0.8,0.8,[],options);
        k = k + 1;
    end
    thetanew = interp1(1:Dt:T,thetaEM,1:T,'linear',thetaEM(end)); % on interpole sur tous les echantillons
    
    nllold = nll;
    nll = negloglikelihoodsig(y,Sigmay); % neg log-vraisemblance qui decroit
    nllV(nbit) = nll;
    
    stopcrit = nllold - nll;
    
    if length(varargin) == 1
        errEM = sum( abs( thetanew(:) - theta(:) ).^2 );
        fprintf(' Iteration %i \n Quadratic Error: %.3f \n\n', nbit, errEM)
    else
        fprintf(' Iteration %i \n Negative loglikelihood: %.3f \n\n', nbit, nllV(nbit))
    end
    
    %plot(thetanew); drawnow;
    nbit = nbit + 1;
    thetaold = thetanew;
end
close; 

thetaEST = interp1(1:Dt:T,thetaEM,1:T,'linear',thetaEM(end));
dgammaEST = 2.^thetaEST;
x = statAMWP(y,ones(T,1),dgammaEST);
SxEST = estim_spec(x,T,alpha); % spectre
W = transform_adap(y,sigmay,TT,Delta,M_psi,Sxest,thetaold,MatPsi);