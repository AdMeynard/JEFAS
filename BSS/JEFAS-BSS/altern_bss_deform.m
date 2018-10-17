function [Boptim,dgamma,Sx,errSDR,errSIR,errSAR,estSIR] = altern_bss_deform(yv,z,B0,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramS,stop_crit,cvmax,stopSIR,options)

[N,T] = size(z);

wav_typ = cell2mat(paramWAV(1));
wav_param = cell2mat(paramWAV(2));

Nf = cell2mat(paramS(2));

scalesBSS = cell2mat(paramBSS(1));
tau0 = cell2mat(paramBSS(2));
dtau = cell2mat(paramBSS(3));
r = cell2mat(paramBSS(4));

Ms = length(scalesBSS);
% Wavelet transforms of the observations:
for n=1:N
    Wz(:,:,n) = cwt(z(n,:),scalesBSS,wav_typ,wav_param);
end

M_psi = bas_calc_cov(scalesBSS,wav_typ,wav_param,2*Nf-1);

% Initializations:
Boptim = B0;
dgamma = dgamma0;
C = zeros(Ms,Ms,N);

nonlcon = @norm_row; % Bien ?

hatyold = z;
% Bold = zeros(N);
Bold = eye(N);

haty = Boptim*z;
haty = reordersig(hatyold,haty);
[SDR,SIR,SAR,~]=bss_eval_sources(haty,yv);
cv = 1;
SIRit = 0;
while ((cv<=cvmax)&&(mean(SIRit)<=stopSIR))
    
    % Warping and spectrum estimation
    dgammaprec = dgamma;
    for n=1:N
        y = haty(n,:);
        [~,dgammaML, Sxn] = estim_altern(y,Dt,dgammaprec(n,:),ones(1,T),paramWAV,paramWP,{'no AM'},paramS,stop_crit,1);
        dgamma(n,:) = dgammaML;
        Sx(n,:) = Sxn;
    end
    
    % Mixing matrix estimation
%     l = 1;
%     for ttau = tau0:(tau0+dtau)
%         for n=1:N
%             Cn = calc_cov(M_psi,Sx(n,:),log2(dgamma(n,ttau)));
%             C(:,:,n,l) = (1-r)*Cn + r*eye(Ms); % regularisation
%         end
%         l = l+1;
%     end
%     costB = @(B)cost_lh2(B,tau0,dtau,Wz,C);
    for n=1:N
        Cn = calc_cov(M_psi,Sx(n,:),log2(dgamma(n,tau0)));
        C(:,:,n) = (1-r)*Cn + r*eye(Ms); % regularisation
    end
    costB = @(B)cost_lh3(B,tau0,Wz,C);
    Boptim = fmincon(costB,Bold,[],[],[],[],-ones(N),+ones(N),nonlcon,options);
    Bold = Boptim;
    
    % reordering sources to keep correspondence between dgamma and haty
    oldhaty = haty;
    haty = Boptim*z;
    haty = reordersig(oldhaty,haty);

    % convergence ?
    SDR_old = SDR;
    SIR_old = SIR;
    SAR_old = SAR;
    [~,SIRit,~,~]=bss_eval_sources(oldhaty,haty); % Perf en supposant qu'on a trouvé les vraies sources
    [SDR,SIR,SAR,~]=bss_eval_sources(haty,yv); % perfs idéales mais il faut connaitre les vraies sources !!!
    errSDR(cv) = sum(SDR - SDR_old)/N;
    errSIR(cv) = sum(SIR - SIR_old)/N;
    errSAR(cv) = sum(SAR - SAR_old)/N;
    estSIR(cv) = sum(SIRit)/N;
    fprintf(' Iteration %i \n update SDR: %.2f dB \n update SIR: %.2f dB\n update SAR: %.2f dB\n Stopping criterion SIR: %.2f dB\n\n',[cv errSDR(cv) errSIR(cv) errSAR(cv) estSIR(cv)]);
    
    cv = cv+1;
end
