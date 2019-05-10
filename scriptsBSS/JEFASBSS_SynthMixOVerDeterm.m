clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-BSS'));

load('../signals/synthetic_NonstatMixtureCrossingOverDeterm.mat');
% load('../signals/synthetic_NonstatMixtureCrossingSameSpec.mat');
Fs = 44100;

[No,T] = size(z);

Delta = 0.01*T;
thres = 0.5;
N = estimSourcesNb(z,Delta,thres);

%% JEFAS-BSS estimation

dgamma0 = ones(N,T);

Dt = 200;
Dtp = 4000; % pSOBI

Kmat = 200; % number of instants where we estimate the unmixing matrix
vectau = floor(linspace(1,T-1,Kmat)); % corresponding instants
L_bss = 10; 

wav_typ = 'sharp';
wav_param = 500;
wav_paramWP = 20;

NbScales = 100;
scales = 2.^(linspace(0,5,NbScales));
subrateWP = 8; % subsampling step for the scales to ensure the covariance invertibility
scalesWP = scales(1:subrateWP:end);
subrateBSS = 2;
scalesBSS = scales(1:subrateBSS:end);

rAM = 1e-3;

NbScalesS = 110;
scalesS = 2.^(linspace(-1,7,NbScalesS));

paramBSS = {scalesBSS,vectau,L_bss};
paramWAV = {wav_typ,wav_param,wav_paramWP};
paramWP = {scalesWP,0.05};
paramAM = {'no AM'};
% paramAM = {'AM',scales,rAM};
paramS = {scalesS};

stop_crit = 5e-3;
stopSIR = 75; % stopping criterion

combinE = nchoosek(1:No,N);
nbBSS = size(combinE,1);
init_meth = 'sobi';
for indBSS = 1:nbBSS
    combin = combinE(indBSS,:);
    zcom = z(combin,:);
    [heapB_init, vectau_init] = JEFASBSSinit(zcom, init_meth, Dtp);

    tic;
    [heapBoptim,dgammaML,aML,SxML,convergeSIR] = altern_bss_deform_nonstat(zcom,heapB_init,vectau_init,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,stopSIR);
    toc;
    hB{indBSS} = heapBoptim;
    dg{indBSS} = dgammaML;
    SxE{indBSS} = SxML;
    cS{indBSS} = convergeSIR;
    
    hy{indBSS} = nonstatunmixing(zcom,heapBoptim,vectau); % unmixing
end

for indBSS = 1:nbBSS
    [hy{indBSS}, st{indBSS}] = reordersig(y,hy{indBSS});
end

%% Analysis
close all;
t = linspace(0,(T-1)/Fs,T);
nu0 = Fs/4;
freqdisp = [8 4 2 1 0.5]; % Displayed frequencies in kHz
sdisp = log2(nu0./(1e3*freqdisp)); % corresponding log-scales

figure; title('Sources and observations');
for n=1:No
    if n <= N
        Wy = cwt_JEFAS(y(n,:),scales,wav_typ,wav_param);
        subplot(2,No,n); imagesc(t,log2(scales),abs(Wy)); yticks(sdisp); yticklabels(freqdisp);
        ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
    end
    Wz = cwt_JEFAS(z(n,:),scales,wav_typ,wav_param);
    subplot(2,nbBSS,No+n); imagesc(t,log2(scales),abs(Wz)); yticks(sdisp); yticklabels(freqdisp);
    xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
end

figure; title('Estimated sources');
for k = 1:nbBSS
    haty = hy{k};
    for n=1:N
        Why = cwt_JEFAS(haty(n,:),scales,wav_typ,wav_param);
        subplot(nbBSS,N,(k-1)*N+n); imagesc(t,log2(scales),abs(Why));
        yticks(sdisp); yticklabels(freqdisp);
        xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
    end
end

figure; title('Estimated Time warpings');
for k = 1:nbBSS
    dgammaML = dg{k};
    stablematch = st{k};
    for n=1:N
        subplot(nbBSS,N,(k-1)*N+n); plot(t,dgamma(n,:),'b--',t,dgammaML(stablematch(n),:),'r','linewidth',2);
        xlabel('Time (s)'); ylabel('\gamma''(t)'); grid on; legend('Ground truth function','Estimated function'); set(gca,'fontsize',18);
    end
end

figure; title('Estimated spectra');
for k = 1:nbBSS
    haty = hy{k};
    dgammaML = dg{k};
    stablematch = st{k};
    for n=1:N
        hatx = statAMWP(haty(n,:),aML(stablematch(n),:),dgammaML(stablematch(n),:));
        alpha = 15;
        Nff = 50000;
        Sxw = estim_spec(hatx,Nff,alpha);
        Sxw = Sxw*sum(Sx(n,:)  + 25e-4)/sum(Sxw); % normalization
        freq = linspace(0,Fs,Nff);
        freq2 = linspace(0,(T-1)*Fs/T,Fs);
        subplot(nbBSS,N,(k-1)*N+n); plot(freq2/1e3,log10(Sx(n,:)  + 25e-4),'b--',freq/1e3,log10(Sxw),'r','linewidth',2); 
        xlabel('Frequency (kHz)'); ylabel('S_x'); axis tight; grid on;
        xlim([0 13]); ylim([-5.1 1]); yticks([-5 -4 -3 -2 -1 0]); yticklabels({'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '10^{0}'}); 
        legend('Ground truth spectrum','Estimated spectrum'); set(gca,'fontsize',18);
    end
end

%% Performances
heapBtot = zeros(N,No,Kmat);

figure;
for indBSS = 1:nbBSS
    combin = combinE(indBSS,:);
    zcom = z(combin,:);
    
    heapBoptim = hB{indBSS};
    haty = hy{indBSS}; % unmixing
    
    [SDR{indBSS},SIR{indBSS}] = bss_eval_sources(haty,y); % JEFAS-BSS
    for k=1:Kmat
       indJEFAS(k) = amari(heapBoptim(:,:,k),heapA(combin,:,vectau(k)));
       heapBtot(:,combin,k) = heapBtot(:,combin,k) + heapBoptim(st{indBSS},:,k);
    end
    indam{indBSS} = indJEFAS;
    
    fprintf('Index   | JEFAS-BSS \n')
    fprintf('SIR     |  %.2f\n', mean(SIR{indBSS}))
    fprintf('SDR     |  %.2f\n', mean(SDR{indBSS}))
    fprintf('Amari   |  %.2f\n', mean(indJEFAS))
    
    hold on; plot(t(vectau),indJEFAS,'linewidth',2); grid on; axis tight; 
    xlabel('Time (s)'); ylabel('Amari index (dB)'); grid on; set(gca,'fontsize',24);
end

%% BSS moyenne
heapBtot = heapBtot/nchoosek(No-1,N-1);
hytot = nonstatunmixing(z,heapBtot,vectau); % unmixing

[SDRtot,SIRtot] = bss_eval_sources(hytot,y); % JEFAS-BSS

for k=1:Kmat
   indJEFAStot(k) = amari(heapBtot(:,:,k),heapA(:,:,vectau(k)));
end

fprintf('Index   | JEFAS-BSS \n')
fprintf('SIR     |  %.2f\n', mean(SIRtot))
fprintf('SDR     |  %.2f\n', mean(SDRtot))
fprintf('Amari   |  %.2f\n', mean(indJEFAStot))


