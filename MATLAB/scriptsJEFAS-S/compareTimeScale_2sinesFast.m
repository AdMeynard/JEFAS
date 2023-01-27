%% Comparison between time-scale representations of the narrowband signal

%% Load signal and JEFAS-S results

clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath('../JEFAS-S/');

load('../signals/TwoSineWaves_fast.mat');
T = length(y);
Fs = T ;
t = linspace(0,(T-1)/Fs,T) ;

load('results/JEFASS_2sinesfast') ;

%% Show CWT

Ms = 1000 ;
fmin = floor(0.048*T) ;
fmax = floor(0.40*T) ;

wav_typ = 'sharp'; % wavelet type (cf. cwt_JEFAS.m)
wav_param = 50 ;
freqtick = [50 100 200 300 400];

figure;

subplot(5,2,[1 9]); 
[Wy,t,s] = display_cwt_JEFAS(y,Fs,fmin,fmax,Ms,wav_typ,wav_param,freqtick) ;
colormap(1-gray) ;
ylabel('Frequency (Hz)');
set(gca,'fontsize',18);

%% Show ordinary Synthesis Time-Scale Representation

prior = 'wavelet' ;
scales = 2.^s ;

W = SynthTSrep(y,sigmay,prior,scales,wav_typ,wav_param,SxEST,log2(dgammaEST)) ;

freqtick = sort(freqtick,'descend');
xi0 = Fs/4; % wavelet central frequency
sdisp = log2(xi0./freqtick); % corresponding log-scales

subplot(5,2,[2 10]); 
imagesc(t,log2(scales),log1p(abs(W)/0.01));
yticks(sdisp); yticklabels(freqtick);
colormap(1-gray) ;
ylabel('Frequency (Hz)');
set(gca,'fontsize',18);

% Show sharp Synthesis Time-Scale Representation

thres = 0.05 ;
[SxSp,freqMax] = sparsifySpectrum(SxEST,thres) ;

fest = freqMax * Fs ;

prior = 'wavelet' ;
wav_param = 1000 ;

Wsharp = SynthTSrep(y,sigmay,prior,scales,wav_typ,wav_param,SxSp,log2(dgammaEST));

figure;

imagesc(t,log2(scales),log1p(abs(Wsharp)/0.1));
yticks(sdisp); yticklabels(freqtick);
colormap(1-gray) ;
xlabel('Time (s)') ;
ylabel('Frequency (Hz)');
set(gca,'fontsize',18);

%% Synchrosqueezing and Second-Order Synchrosqueezing

addpath('SST2/')

No = floor(log2(T))-1;
Nvpo = min(floor(round(length(scales)/No) / 2)*2,48) ;

gamma = 0.01 ;
mywav = 'cmor2-1' ;
[Wav,WavSST,WavSST2,f,sW, omega, omega2, tau, phipp, norm2psi] = Wsst2(y,gamma,mywav,Nvpo) ;

figure;
imagesc(t,log(f),log1p(abs(WavSST)/0.1)); axis xy;
yticks(sort(log(freqtick))); yticklabels(sort(freqtick)) ;
ylim([log(fmin) log(fmax)])
colormap(1-gray) ;
xlabel('Time (s)') ;
ylabel('Frequency (Hz)');
set(gca,'fontsize',18);

figure;
imagesc(t,log(f),log1p(abs(WavSST2)/0.1)); axis xy;
yticks(sort(log(freqtick))); yticklabels(sort(freqtick)) ;
ylim([log(fmin) log(fmax)])
colormap(1-gray) ;
xlabel('Time (s)') ;
ylabel('Frequency (Hz)');
set(gca,'fontsize',18);

% %% Adaptive WSST
% 
% addpath('Toolbox_ASST_V1/');
% 
% mu = 1 ;
% sgm = linspace(1,5,200) ; % potential values for sigma
% t1 = 16 ; % number of points for integration
% sig_u = Renyi_entropy(y,sgm,2.5,t1,mu,Nvpo) ;
% 
% [Ws, Ts1_vt, Ts2_vt, w2nd, aj] = cwt_2nd_order_sst_tv(y,gamma,mu,sig_u,Nvpo);
% figure;
% imagesc(abs(Ts1_vt));
% figure;
% imagesc(abs(Ts2_vt));
