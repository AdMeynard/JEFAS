%% JEFAS-S initialized with the results of JEFAS

clear all; close all; clc;
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-S'));
addpath('../signals');
addpath('../cwt');

load('toysigJEFASS')
theta = log2(dgamma);
T = length(y0);
%% JEFAS-S

Nbscales = 100;
scales = 2.^linspace(0,3,Nbscales);
scales = scales(1:7:end); % pour se comparer a JEFAS
wav_param = 300;
wav_paramWP = 10; % used by JEFAS and JEFASS
wav_typ = 'sharp';

% parameters specified for initialization by JEFAS algorithm
ratio = 0.05; % activity threshold
dgamma0 = ones(1,T);
a0 = ones(1,T);
stop_crit = 5e-3; % relative update threshold
paramWAV = {wav_typ,wav_param,wav_paramWP};
paramWP = {scales};
paramAM = {'no AM'};
scalesS = 2.^(linspace(-1,3.5,150)); % for spectrum estimation
paramS = {scalesS};

itD = 10; % nb d'iterations de la descente de gradient
stopD = 1e-4; % tolerance de la descente de gradient
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxIterations',itD,'StepTolerance',stopD,'Display','off'); % fmincon

Dt = 50; % valeur du sous echantillonage pour l'estimation de thetaEM
TT = 512; %1024;%T/16;%256 % taille de la decoupe pour W
Delta = 0.75*TT; % recouvrement
alpha = 15; % largeur pour l'estimation du spectre par Welch

Nit = 5; % number of iterations
thres = 0.2; % ?

wav_paramR = 200;
scalesR = 2.^linspace(0,3,Nbscales);
[M_psi,M_tmpdpsi] = bas_calc_dcov(scalesR,wav_typ,wav_paramR,T);
MatPsi = ifft(M_psi.',[],1);

sigmaH = [5e-5 5e-4 5e-3 5e-2 5e-1 1];
P = length(sigmaH);

for p = 1:P
    
    sigmay = sigmaH(p) ;
    y = y0 + sigmay*randn(T,1) ;
    
    %%% APPLY JEFAS INITIALISATION
    [~, dgammaML, ~, ~] = estim_altern(y,Dt,ratio,dgamma0,a0,paramWAV,paramWP,paramAM,paramS,stop_crit,10,0);
    thetaINIT = log2(dgammaML); % Initialization from JEFAS results
    
    % Apply JEFAS-S
    tic;
    [dgammaEST, SxEST, W, nll] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_paramWP,Dt,TT,Delta,alpha,Nit,thres,itD,stopD,theta);
    toc;

    dgammaH{p} = dgammaEST ;
    SxESTH{p} = SxEST ;
    WH{p} = W ;

    uEST = abs(dgammaEST(:)-dgamma).^2;
    uML = abs(dgammaML(:)-dgamma).^2;

    fprintf('Noise variance %.2e\n\n', sigmaH(p)^2);
    fprintf('EQM:\n JEFAS: %.3e\n JEFAS-S: %.3e\n\n',mean(uML(1:8151)),mean(uEST(1:8151)));  % prevent edge effect
    
    [Wfinest, MMSigmay] = transform_adap(y,sigmay,TT,Delta,M_psi,SxEST,log2(dgammaEST),MatPsi); % adpated representation
    Cy = buildSigmay(MMSigmay,T,TT,Delta);
    yr = synthesis(Wfinest, MatPsi); % reconstruction

    snr0th = 10*log10(1 + var(y0)/sigmay^2);
    snr0 = snr(y,y-y0);
    snrc = snr(yr(:),yr(:)-y0);
    fprintf('SNR original %.2f\n SNR wavelet-like %.2f\n\n',snr0,snrc)
    
    biasTH = theo_bias(y0,TT,Delta,MMSigmay,sigmay);
    biasXP = yr(:)-y0;

    Cyi = inv(Cy);
    varty0 = trace( (eye(T)-sigmay^2*Cyi)^2 * (sigmay^2*eye(T) + y0*y0.') );
    vartnoise = sigmay^2 * ( trace((eye(T)-sigmay^2*Cyi)^2) + trace(sigmay^2*Cyi^2*y0*y0.') );
    SNR_TH = 10*log10( varty0 / vartnoise ); % Theoretical SNR
    fprintf('Theoretical SNR improvment %.2f\n Experimental SNR improvment %.2f\n\n',SNR_TH-snr0th,snrc-snr0)

    fprintf('------------------------------------------------\n\n')
end

%% Results and Adapted representation

% t = linspace(0,(T-1)/Fs,T);
% figure;
% plot(nll,'linewidth',2); grid on;
% xlabel('Iteration'); ylabel('negative log-likelihood of the signal');
% 
% figure;
% plot(t,dgamma,'b--',t,dgammaML,'k',t,dgammaEST,'r--','linewidth',2); grid on;
% xlabel('Time (s)'); ylabel('\gamma''(t)')
% legend({'Ground truth fonction','JEFAS estimate','JEFAS-S estimate'},'FontSize',20);
% set(gca,'FontSize',20);

% Wfin = transform_adap(y,sigmay,TT,Delta,M_psi,Sx,log2(dgamma),MatPsi);

% xi0 = Fs/4; % wavelet central frequency
% freqdisp = [2 1 0.5 0.3]; % Displayed frequencies in kHz
% sdisp = log2(xi0./(freqdisp*1e3));
% 
% Wy = cwt_JEFAS(y,scales,wav_typ,wav_param);
% 
% figure;
% subplot(2,1,1);
% plot(t,y,'linewidth',2);
% xlabel('Time (s)'); ylabel('Signal');
% set(gca,'FontSize',20);
% subplot(2,1,2);
% imagesc(t,log2(scales),abs(Wy));
% xlabel('Time (s)'); ylabel('Frequency (kHz)');
% %colormap(1-gray); 
% yticks(sdisp); yticklabels(freqdisp);set(gca,'FontSize',20);
% 
% figure;
% subplot(1,2,2);
% imagesc(t,log2(scales),abs(Wfinest));
% xlabel('Time (s)'); ylabel('Frequency (kHz)');
% yticks(sdisp); yticklabels(freqdisp);
% set(gca,'FontSize',20);
% subplot(1,2,1);
% imagesc(t,log2(scales),abs(Wy));
% xlabel('Time (s)'); ylabel('Frequency (kHz)');
% %colormap(1-gray); 
% yticks(sdisp); yticklabels(freqdisp);set(gca,'FontSize',20);

%% Synthesis

% figure; plot(t,y0,'b--',t,y,'k',t,yr,'r','linewidth',2)

% figure;
% plot(t, biasTH, t, biasXP, 'linewidth', 2); axis tight; 
% xlim([0 1]); grid on; 
% legend({'Theoretical bias', 'Experimental bias'},'FontSize',20);
% set(gca,'FontSize',20);