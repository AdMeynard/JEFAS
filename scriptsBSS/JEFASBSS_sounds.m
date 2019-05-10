clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-BSS'));

load('../signals/sound_Mixture.mat');
[N,T] = size(z);

%% SOBI estimation
p = 1000; % Number of jointly diagonalized matrices
[Asobi,ysobi] = sobi(z,N,p);
Bsobi = inv(Asobi);

%% p-SOBI estimation
Dtp = 4000;
vectauns = 1:Dtp:T;
[heap_Apsobi, heap_Bpsobi] = psobi(z,vectauns);
haty0 = nonstatunmixing(z,heap_Bpsobi,vectauns);

%% JEFAS-BSS estimation

dgamma0 = ones(N,T);

Dt = 500;

Kmat = 300; % number of instants where we estimate the unmixing matrix
vectau = floor(linspace(1,T-1,Kmat)); % corresponding instants
eps_bss = 10; %0.01; %

wav_typ = 'sharp';
wav_param = 1000;
wav_paramWP = 20;

NbScales = 125;
scales = 2.^(linspace(1,5.5,NbScales));
subrateWP = 3; % subsampling step for the scales to ensure the covariance invertibility
scalesWP = scales(1:subrateWP:end);
subrateBSS = 2;
scalesBSS = scales(1:subrateBSS:end);

rAM = 1e-4;

L_bss = 10; % nb it newton

NbScalesS = 110;
scalesS = 2.^(linspace(-1,7,NbScalesS));

paramBSS = {scalesBSS,vectau,L_bss};
paramWAV = {wav_typ,wav_param,wav_paramWP};
paramWP = {scalesWP,0.1};
paramAM = {'no AM'};
% paramAM = {'AM',scales,rAM};
paramS = {scalesS};

stop_crit = 5e-3;
stopSIR = 75; % stopping criterion

Nit = 7;

init_meth = 'sobi';
[heapB_init, vectau_init] = JEFASBSSinit(z, init_meth, Dtp);

tic;
[heapBoptim,dgammaML,aML,SxML,convergeSIR] = altern_bss_deform_nonstat(z,heapB_init,vectau_init,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,stopSIR,Nit);
toc;

haty = nonstatunmixing(z,heapBoptim,vectau); % unmixing


%% Analysis
close all;
t = linspace(0,(T-1)/Fs,T);
nu0 = Fs/4;
freqdisp = [16 8 4 2 1 0.5 0.25]; % Displayed frequencies in kHz
sdisp = log2(nu0./(1e3*freqdisp)); % corresponding log-scales

figure; title('Sources and observations');
for n=1:N
    Wy = cwt_JEFAS(y(n,:),scales,wav_typ,wav_param);
    subplot(N,N,n); imagesc(t,log2(scales),abs(Wy)); yticks(sdisp); yticklabels(freqdisp);
    ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
    Wz = cwt_JEFAS(z(n,:),scales,wav_typ,wav_param);
    subplot(N,N,N+n); imagesc(t,log2(scales),abs(Wz)); yticks(sdisp); yticklabels(freqdisp);
    xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
end

[haty, stablematch] = reordersig(y,haty);

figure; title('Estimated sources');
for n=1:N
    Why = cwt_JEFAS(haty(n,:),scales,wav_typ,wav_param);
    subplot(1,N,n); imagesc(t,log2(scales),abs(Why));
    yticks(sdisp); yticklabels(freqdisp);
    xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
end

figure; title('Estimated Time warpings');
for n=1:N
    hold on; plot(t,dgammaML(stablematch(n),:),'linewidth',2);
    xlabel('Time (s)'); ylabel('\gamma''(t)'); grid on; legend('Estimated wolf 1','Estimated wolf 2'); set(gca,'fontsize',18);
end

figure; title('Estimated spectra');
for n=1:N
    hatx = statAMWP(haty(n,:),aML(stablematch(n),:),dgammaML(stablematch(n),:));
    alpha = 15;
    Nff = 50000;
    Sxw = estim_spec(hatx,Nff,alpha);
    freq = linspace(0,Fs,Nff);
    freq2 = linspace(0,(T-1)*Fs/T,Fs);
    hold on; plot(freq/1e3,log10(Sxw),'linewidth',2); 
    xlabel('Frequency (kHz)'); ylabel('Sprectra'); axis tight; grid on;
    xlim([0 5]); 
%     ylim([-5.1 1]); yticks([-5 -4 -3 -2 -1 0]); yticklabels({'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '10^{0}'}); 
    legend('Source 1','Source 2'); set(gca,'fontsize',18);
end

%% Performances
[SDRsobi,SIRsobi] = bss_eval_sources(ysobi,y); % SOBI
[SDRpsobi,SIRpsobi] = bss_eval_sources(haty0,y); %p-SOBI
% [SDRqtf,SIRqtf] = bss_eval_sources(hatyqtf,y); % QTF-BSS
[SDR,SIR] = bss_eval_sources(haty,y); % JEFAS-BSS

% Amari index
for k=1:Kmat
   indJEFAS(k) = amari(heapBoptim(:,:,k),heapA(:,:,vectau(k)));
   indSOBI(k) = amari(Bsobi,heapA(:,:,vectau(k)));
%    indQTF(k) = amari(BQTF,heapA(:,:,vectau(k)));
   k2 = length(vectauns(vectauns<=vectau(k)));
   indPSOBI(k) = amari(heap_Bpsobi(:,:,k2),heapA(:,:,vectau(k)));
end

fprintf('Index   |   SOBI  |  p-SOBI | JEFAS-BSS \n')
fprintf('SIR     |  %.2f  |   %.2f  |  %.2f\n', mean(SIRsobi),mean(SIRpsobi),mean(SIR))
fprintf('SDR     |  %.2f  | %.2f  |  %.2f\n', mean(SDRsobi),mean(SDRpsobi),mean(SDR))
fprintf('Amari   |  %.2f  |  %.2f  | %.2f\n', mean(indSOBI),mean(indPSOBI),mean(indJEFAS))

figure;subplot(2,1,2);
plot(t(vectau),indSOBI,'k--',t(vectau),indPSOBI,'r:',t(vectau),indJEFAS,'b','linewidth',2); grid on; axis tight; 
xlabel('Time (s)'); ylabel('Amari index (dB)'); grid on; legend('SOBI','p-SOBI','JEFAS-BSS'); set(gca,'fontsize',24);
%% Convergence

subplot(2,1,1);
plot(convergeSIR,'linewidth',2);
xlabel('Iteration'); ylabel('Convergence criterion (dB)'); grid on; axis tight;