%% Algo EM en initialisant par le resultat de JEFAS

clear all; close all; clc;
addpath(genpath('../JEFAS-S'));
addpath('../signals');
addpath('data');

% %% Synthèse d'un signal large bande loc. dilaté
% Fs = 2^13; %Fs = 2^10;
% T = Fs;
% t = linspace(0,(T-1)/Fs,T);
% 
% x = randn(T,1);
% [x1,Sx1] = BandPassApprox(x,floor(0.10*T),floor(0.12*T));
% [x2,Sx2] = BandPassApprox(x,floor(0.05*T),floor(0.08*T));
% x = x1+x2;
% Sx = Sx1 + Sx2;
% 
% % [y0,~,dgamma] = sinewarp(x,2,0.02);
% [y0,gamma,dgamma] = newarp(x,2,2);
% 
% sigmay = 5e-2;
% y = y0 + sigmay*randn(T,1);

load('toysig_short')
theta = log2(dgamma);
T = length(y);
%% algo EM

Nbscales = 100;
scales = 2.^linspace(0,3,Nbscales);
scales = scales(1:7:end); % pour se comparer a JEFAS
wav_param = 10;
wav_typ = 'sharp';

itD = 10; % nb d'iterations de la descente de gradient
stopD = 1e-4; % tolerance de la descente de gradient
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxIterations',itD,'StepTolerance',stopD,'Display','off'); % fmincon

Dt = 50; % valeur du sous echantillonage pour l'estimation de thetaEM
TT = 512; %1024;%T/16;%256 % taille de la decoupe pour W
Delta = 0.75*TT; % recouvrement
alpha = 15; % largeur pour l'estimation du spectre par Welch

Nit = 5; % nombre d'iterations de EM

initMeth = 'JEFAS';
switch initMeth
    
    case 'true'
        thetaINIT = log2(dgamma);
    
    case 'none'
        thetaINIT = zeros(1,T);
        
    case 'JEFAS'
        load('resJEFAS_toysig.mat')
        thetaINIT = log2(dgammaML); % Initialisation a partir des resultats de JEFAS
        
    case 'gwn'
        [M_psi,M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,T);
        MatPsi = ifft(M_psi.',[],1);
        [W, MMSigmay] = transform_adap(y,sigmay,TT,Delta,M_psi,repmat(var(y),1,T),zeros(1,T),MatPsi);
        Sigmay = buildSigmay(MMSigmay,T,TT,Delta); % matrice de tout le signal
        iSigmay = inv(Sigmay);
        Syest = estim_spec(y,T,alpha);
        k = 1;
        for t = 1:Dt:T
            U = W(:,t);
            theta0 = 0;
            Qemx = @(x)Qem(x, theta0, t, U, M_psi, M_tmpdpsi, Syest, MatPsi, iSigmay); % fonction a minimiser
            thetaEM(k) = fmincon(Qemx,theta0,[],[],[],[],-0.8,0.8,[],options);
            k = k + 1;
        end
        thetaINIT = interp1(1:Dt:T,thetaEM,1:T,'linear',thetaEM(end)); % on interpole sur tous les echantillons
end

thres = 0.2;
tic;
[dgammaEST, SxEST, W, nll] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_param,Dt,TT,Delta,alpha,Nit,thres,itD,stopD,theta);
toc;

%% Results
t = linspace(0,(T-1)/Fs,T);

figure;
plot(nll,'linewidth',2); grid on;
xlabel('Iteration'); ylabel('neg log-vraisemblance du signal');

figure;
plot(t,dgamma,'b--',t,dgammaML,'k',t,dgammaEST,'r--','linewidth',2); grid on;
xlabel('Time (s)'); ylabel('\gamma''(t)')
legend({'Ground truth fonction','JEFAS estimate','JEFAS-S estimate'},'FontSize',20);
set(gca,'FontSize',20);

uEST = abs(dgammaEST(:)-dgamma).^2;
uML = abs(dgammaML(:)-dgamma).^2;
fprintf('EQM:\n JEFAS: %.3e\n JEFAS-S: %.3e\n',mean(uML(1:8151)),mean(uEST(1:8151))) % prevent edge effect
%% Adapted representation
wav_param = 200;
scales = 2.^linspace(0,3,Nbscales);
[M_psi,M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,T);
MatPsi = ifft(M_psi.',[],1);

% Wfin = transform_adap(y,sigmay,TT,Delta,M_psi,Sx,log2(dgamma),MatPsi);
[Wfinest, MMSigmay] = transform_adap(y,sigmay,TT,Delta,M_psi,SxEST,log2(dgammaEST),MatPsi);

xi0 = Fs/4; % wavelet central frequency
freqdisp = [2 1 0.5 0.3]; % Displayed frequencies in kHz
sdisp = log2(xi0./(freqdisp*1e3));

Wy = cwt_JEFAS(y,scales,wav_typ,wav_param);

figure;
subplot(2,1,1);
plot(t,y,'linewidth',2);
xlabel('Time (s)'); ylabel('Signal');
set(gca,'FontSize',20);
subplot(2,1,2);
imagesc(t,log2(scales),abs(Wy));
xlabel('Time (s)'); ylabel('Frequency (kHz)');
%colormap(1-gray); 
yticks(sdisp); yticklabels(freqdisp);set(gca,'FontSize',20);

figure;
subplot(1,2,2);
imagesc(t,log2(scales),abs(Wfinest));
xlabel('Time (s)'); ylabel('Frequency (kHz)');
yticks(sdisp); yticklabels(freqdisp);
set(gca,'FontSize',20);
subplot(1,2,1);
imagesc(t,log2(scales),abs(Wy));
xlabel('Time (s)'); ylabel('Frequency (kHz)');
%colormap(1-gray); 
yticks(sdisp); yticklabels(freqdisp);set(gca,'FontSize',20);

Cy = buildSigmay(MMSigmay,T,TT,Delta);
figure;
imagesc(t,t,abs(Cy));
xlabel('Time (s)'); ylabel('Time (s)');
set(gca,'FontSize',20); %colormap(1-gray); 
title('Signal expected covariance matrice');

%% Synthesis
yr = synthese(Wfinest, MatPsi);
figure; plot(t,y0,'b--',t,y,'k',t,yr,'r','linewidth',2)

snr0th = 10*log10(1 + var(y0)/sigmay^2);
snr0 = snr(y,y-y0);
snrc = snr(yr(:),yr(:)-y0);
fprintf('SNR original %.2f\n SNR wavelet-like %.2f\n',snr0,snrc)

biasTH = theo_bias(y0,TT,Delta,MMSigmay,sigmay);
biasXP = yr(:)-y0;

figure;
plot(t, biasTH, t, biasXP, 'linewidth', 2); axis tight; 
xlim([0 1]); grid on; 
legend({'Theoretical bias', 'Experimental bias'},'FontSize',20);
set(gca,'FontSize',20);

Cyi = inv(Cy);
varty0 = trace( (eye(T)-sigmay^2*Cyi)^2 * (sigmay^2*eye(T) + y0*y0.') );
vartnoise = sigmay^2 * ( trace((eye(T)-sigmay^2*Cyi)^2) + trace(sigmay^2*Cyi^2*y0*y0.') );
SNR_TH = 10*log10( varty0 / vartnoise ); % Theoretical SNR
fprintf('Theoretical SNR improvment %.2f\n Experimental SNR improvment %.2f\n',SNR_TH-snr0th,snrc-snr0)
