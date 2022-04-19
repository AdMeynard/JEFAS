clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath('../JEFAS-S/');

%% Load signal and JEFAS-S results

load('results/JEFAS_wind4JEFASS') ;

load('../signals/wind4JEFASS') ;
T = length(y);
Fs = T ;
t = linspace(0,(T-1)/Fs,T) ;
f = linspace(0,(T-1)*Fs/T,T);

load('results/JEFASS_wind4JEFASS') ;

%% Show JEFAS-S output

figure;
plot(t,2.^dgamma,'k',t,2.^dgammaML,'b',t,2.^dgammaEST,'r','linewidth',2);
xlabel('Time (s)'); ylabel("\gamma'")
legend('Ground-truth function','JEFAS estimate','JEFAS-S estimate');
set(gca,'fontsize',18);

figure;
plot(f,Sx,'k',f,SxEST,'r','linewidth',2);
xlim([0 Fs/2]);
xlabel('frequency (Hz)'); ylabel("Power spectrum")
legend('Ground-truth function','JEFAS-S estimate');
set(gca,'fontsize',18);

%% Show ordinary Synthesis Time-Scale Representation

% Scales:
fmin = floor(0.058*T) ;
fmax = floor(0.42*T) ;
xi0 = Fs/4; % wavelet central frequency
smin = log2( xi0/fmax );
smax = log2( xi0/fmin );

Ms = 200 ;
scales = 2.^(linspace(smin,smax,Ms));

wav_typ = 'sharp';
wav_param = 100 ;
prior = 'wavelet' ;

W = SynthTSrep(y,sigmay,prior,scales,wav_typ,wav_param,SxEST,log2(dgammaEST),T) ;

freqtick = [100 200 400 800];
freqtick = sort(freqtick,'descend');
xi0 = Fs/4; % wavelet central frequency
sdisp = log2(xi0./freqtick); % corresponding log-scales
 
figure;
subplot(2,1,1);
imagesc(t,log2(scales),log1p(abs(W)/0.1));
yticks(sdisp); yticklabels(freqtick);
colormap(1-gray) ;
ylabel('Frequency (Hz)');
set(gca,'fontsize',18);

%% Synthesize the signal

[M_psi,~] = bas_calc_dcov(scales,wav_typ,wav_param,T);
MatPsi = ifft(M_psi.',[],1);

yr = synthesis(W, MatPsi); % reconstruction

subplot(2,1,2);
plot(t,y,'k',t,y0,'b',t,yr,'r','linewidth',2); xlim([0 0.2])
xlabel('Time (s)'); ylabel('Signals')
legend('Observed signal','Meaningful signal','Synthesized signal');
set(gca,'fontsize',18);