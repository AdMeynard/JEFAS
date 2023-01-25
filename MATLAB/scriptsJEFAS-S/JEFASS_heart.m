%% Algo EM sur le signal cardiaque

clear all; close all; clc;
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-S'));
addpath(genpath('../cwt'));
addpath('../signals');
%% Signal et resultats de JEFAS
load('SingleComp'); y = s; clear s;
T = length(y);

%% algo EM

Nbscales = 100;
scales = 2.^linspace(-1,3,Nbscales);
scales = scales(1:7:end); % pour se comparer a JEFAS
wav_param = 10;
wav_typ = 'sharp';

prior = 'wavelet' ; %prior for the covariance matrix of the time scale representation
priorList = {prior} ;

itD = 10; % nb d'iterations de la descente de gradient
stopD = 1e-4; % tolerance de la descente de gradient
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxIterations',itD,'StepTolerance',stopD,'Display','off'); % fmincon

Dt = 2; % valeur du sous echantillonage pour l'estimation de thetaEM
TT = 48; % taille de la decoupe pour W
Delta = 0.75*TT; % recouvrement
alpha = 15; % largeur pour l'estimation du spectre par Welch

sigmay = 0.05; % noise variance

TS = T ;

Nit = 75; % nombre d'iterations de EM

initMeth = 'gwn';
switch initMeth
    
    case 'none'
        thetaINIT = ones(1,T);
        
    case 'JEFAS'
        thetaINIT = thetaML; % Initialisation a partir des resultats de JEFAS

        case 'gwn'
        [M_psi,M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,T);
        MatPsi = ifft(M_psi.',[],1);
        [W, MMSigmay] = transform_adap(y,sigmay,priorList,TT,Delta,M_psi,repmat(var(y),1,T),zeros(1,T),MatPsi);
        Sigmay = buildSigmay(MMSigmay,T,TT,Delta); % matrice de tout le signal
        iSigmay = inv(Sigmay);
        Syest = estim_spec(y,T,alpha);
        k = 1;
        for t = 1:Dt:T
            U = W(:,t);
            theta0 = 0;
            Qemx = @(x)Qem(x, theta0, t, U, priorList, M_psi, M_tmpdpsi, Syest, MatPsi, iSigmay); % fonction a minimiser
            thetaEM(k) = fmincon(Qemx,theta0,[],[],[],[],-0.8,0.8,[],options);
            k = k + 1;
        end
        thetaINIT = interp1(1:Dt:T,thetaEM,1:T,'linear',thetaEM(end)); % on interpole sur tous les echantillons

end

thres = 0.4;
tic;
[dgammaEST,SxEST, W, nll] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_param,priorList,Dt,TT,Delta,TS,alpha,Nit,thres,itD,stopD);
toc;

%% Resultats
t = linspace(0,(T-1)/Fs,T);

figure;
plot(nll,'linewidth',2); grid on;
xlabel('Iteration'); ylabel('neg log-vraisemblance du signal');

figure;
plot(t,if1/mean(if1),'b--',t,dgammaEST,'r','linewidth',2); grid on;
xlabel('Time (s)'); ylabel('\gamma''(t)');
legend({'Ground truth fonction','JEFAS-S estimate'},'FontSize',20);
set(gca,'FontSize',20);

%% Adapted representation

% Wavelet-like covariance
wav_param = 15;
scales = 2.^linspace(-1,3,Nbscales);
[M_psi,M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,T);
MatPsi = ifft(M_psi.',[],1);
W = transform_adap(y,sigmay,priorList,TT,Delta,M_psi,SxEST,log2(dgammaEST),MatPsi);

TS = length(SxEST);
omega = (0:(TS-1))*2*pi/TS;
C0 = calc_cov(M_psi,SxEST,0);

% Sparse covariance
sigmaW = 5; %time scale representation concentration around the instantaneous frequency
km = find(SxEST(1:(T/2))==max(SxEST));
S = exp(-((1:(T/2))-km).^2/sigmaW^2);
S = [S fliplr(S)];
S = sum(SxEST)*S/sum(S);
Wnew = transform_adap_sparse(y,sigmay,TT,Delta,scales,S,log2(dgammaEST),MatPsi);

C0new = calc_cov_sparse(scales,S(:),0);

Wy = cwt_JEFAS(y,scales,wav_typ,wav_param); % wavelet transform


xi0 = Fs/4; % wavelet central frequency
freqdisp = [4 2 1 0.5]; % Displayed frequencies in kHz
sdisp = log2(xi0./(freqdisp));

figure;
subplot(1,2,1);
imagesc(log2(scales),log2(scales),abs(C0));
yticks(sdisp); yticklabels(freqdisp);
xticks(sdisp); xticklabels(freqdisp);
xlabel('Frequency (Hz)'); ylabel('Frequency (Hz)');
set(gca,'FontSize',20);
subplot(1,2,2);
imagesc(log2(scales),log2(scales),abs(C0new));
yticks(sdisp); yticklabels(freqdisp);
xticks(sdisp); xticklabels(freqdisp);
xlabel('Frequency (Hz)'); ylabel('Frequency (Hz)');
set(gca,'FontSize',20); %colormap(1-gray);

figure;
subplot(1,2,1);
imagesc(t,log2(scales),log1p(abs(W)));
yticks(sdisp); yticklabels(freqdisp);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'FontSize',20);
subplot(1,2,2);
imagesc(t,log2(scales),log1p(abs(Wnew)));
yticks(sdisp); yticklabels(freqdisp);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'FontSize',20); %colormap(1-gray);

figure;
imagesc(t,log2(scales),log1p(abs(Wy)));
yticks(sdisp); yticklabels(freqdisp);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'FontSize',20); %colormap(1-gray);


yr = synthesis(W, MatPsi);
yrnew = synthesis(Wnew, MatPsi);
figure; plot(t,y,'b--',t,yr,'r',t,yrnew,'k','linewidth',2);
xlabel('Time (s)'); ylabel('Signals');
xlim([10 30])
legend({'Measured signal','Reconstructed signal / covariance wavelet-like','Reconstructed signal / sparse covariance'},'FontSize',20);
set(gca,'FontSize',20);