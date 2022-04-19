%% JEFAS-S initialized with the results of JEFAS

clear all; close all; clc;
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-S'));
addpath('../signals');
addpath('../cwt');
addpath('results');

cd ../JEFAS-S/minFunc_2012       % Change to the unzipped directory
mexAll;               % Compile mex files (not necessary on all systems)

load('toysig_short')
theta = log2(dgamma);
T = length(y);
%% JEFAS-S

Nbscales = 100;
scales = 2.^linspace(0,3,Nbscales);
scales = scales(1:7:end); % for comparison with JEFAS
wav_param = 10;
wav_typ = 'sharp';

prior = 'wavelet' ; %prior for the covariance matrix of the time scale representation
priorList = {prior} ;

itD = 10; % number of iterations for gradient descent
stopD = 1e-4; % gradient descent: tolerance
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxIterations',itD,'StepTolerance',stopD,'Display','off'); % fmincon

Dt = 50; % subsampling value for thetaEM estimation
TT = 512; %  slicing size for W
Delta = 0.75*TT; % overlap
alpha = 15; % width for Welch estimation

TS = 1024 ;

Nit = 5; % number of iterations

initMeth = 'JEFAS';
switch initMeth
    
    case 'true'
        thetaINIT = log2(dgamma);
    
    case 'none'
        thetaINIT = zeros(1,T);
        
    case 'JEFAS'
        load('resJEFAS_toysig.mat')
        thetaINIT = log2(dgammaML); % Initialization from JEFAS results
end

thres = 0.2;
tic;
[dgammaEST, SxEST, W, nll] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_param,priorList,Dt,TT,Delta,TS,alpha,Nit,thres,itD,stopD,theta);
toc;

%% Results
t = linspace(0,(T-1)/Fs,T);

figure;
plot(nll,'linewidth',2); grid on;
xlabel('Iteration'); ylabel('negative log-likelihood of the signal');

figure;
plot(t,dgamma,'b--',t,dgammaML,'k',t,dgammaEST,'r--','linewidth',2); grid on;
xlabel('Time (s)'); ylabel('\gamma''(t)')
legend({'Ground truth fonction','JEFAS estimate','JEFAS-S estimate'},'FontSize',20);
set(gca,'FontSize',20);

uEST = abs(dgammaEST(:)-dgamma).^2;
uML = abs(dgammaML(:)-dgamma).^2;
fprintf('MSE:\n JEFAS: %.3e\n JEFAS-S: %.3e\n',mean(uML(1:8151)),mean(uEST(1:8151))) % prevent edge effect
%% Adapted representation
wav_param = 200;
scales = 2.^linspace(0,3,Nbscales);
[M_psi,M_tmpdpsi] = bas_calc_dcov(scales,wav_typ,wav_param,T);
MatPsi = ifft(M_psi.',[],1);

[Wfinest, MMSigmay] = transform_adap(y,sigmay,priorList,TT,Delta,M_psi,SxEST,log2(dgammaEST),MatPsi);

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
yticks(sdisp); yticklabels(freqdisp);set(gca,'FontSize',20);

Cy = buildSigmay(MMSigmay,T,TT,Delta);
figure;
imagesc(t,t,abs(Cy));
xlabel('Time (s)'); ylabel('Time (s)');
set(gca,'FontSize',20); %colormap(1-gray); 
title('Signal expected covariance matrice');

%% Synthesis
yr = synthesis(Wfinest, MatPsi);
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
