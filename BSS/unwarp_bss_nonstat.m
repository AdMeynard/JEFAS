clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath('JEFAS-BSS');
addpath(genpath('othermethods_BSS')); % quadratic TF BSS mmethod

% load('sig_compBSS.mat'); Fs = 44100; pileA = repmat(A,1,1,Fs);
load('sig_compBSS_nonstat.mat'); Fs = 44100;

[N,T] = size(z);

Kmat = 200; % nombre d'instants ou l'on calcule Boptim
vectau = floor(linspace(1,T-1,Kmat));

%% 1ere methode: via SOBI
p = 1000; % Nombre de matrices à diagonaliser conjointement
[Asobi,ysobi] = sobi(z,N,p);
Bsobi = inv(Asobi);
Bsobi = Bsobi./sqrt(sum(abs(Bsobi).^2,2));

[SDRsobi,SIRsobi,SARsobi,permsobi] = bss_eval_sources(ysobi,y); % SOBI

%% 2eme methode: via SOBI par morceaux
Dt = 4000;
vectauns = 1:Dt:T;
[pile_Asobins, pile_Bsobins] = sobi_nonstat(z,vectauns);
haty0 = nonstatunmixing(z,pile_Bsobins,vectauns);

[SDRsobins,SIRsobins,SARsobins,permsobins] = bss_eval_sources(haty0,y);

%% 3eme methode: via TFQ
eps1 = 10;
eps2 = 0.1;
eps3 = 0.1;
eps4 = 100;

pp = 10; % subsampling 
[~, Dmat] = selecpts(z, pp, eps1, eps2, eps3, eps4); % select_pts
nn = 4;
Aest = BSS_TFQ(Dmat,N,nn); % BSS dessus
BTFQ = inv(Aest);

hatymoreau = BTFQ*z;
[SDRmoreau,SIRmoreau,SARmoreau,permmoreau] = bss_eval_sources(hatymoreau,y);

%% 4eme methode: via mon algo de max de vraisemblance
dgamma0 = dgamma;%ones(N,T);

Dt = 200;

eps_bss = 0.5;
rBSS = 1e-7; %1e-5

wav_typ = 'sharp';
wav_param = 500;
wav_paramWP = 20;

NbScales = 100;
scales = 2.^(linspace(-1,4,NbScales));
subrate = 8; % subsampling step for the scales to ensure the covariance invertibility
scalesWP = scales(1:subrate:end);

rAM = 1e-3;

NbScalesS = 110;
scalesS = 2.^(linspace(-1,7,NbScalesS));

paramBSS = {scales,vectau,eps_bss,rBSS};
paramWAV = {wav_typ,wav_param,wav_paramWP};
paramWP = {scalesWP};
% paramAM = {'no AM'};
paramAM = {'AM',scales,rAM};
paramS = {scalesS};

itMAX = 100;
stopO = 1e-5;
options = optimoptions('fmincon','Algorithm','interior-point','MaxIterations',itMAX,'StepTolerance',stopO,'Display','off');
stop_crit = 5e-3;
cvmax = 20;
stopSIR = 90; % critère d'arret sur l'amélioration

init_meth = 'sobi';
[pileB_init, vectau_init] = JEFASBSSinit(z, init_meth);

tic;
[pileBoptim,dgammaML,aML,Sx,errSIR,errSAR,errSDR,estSIR] = altern_bss_deform_nonstat(y,z,pileB_init,vectau_init,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,cvmax,stopSIR,options);
toc;

haty = nonstatunmixing(z,pileBoptim,vectau); % unmixing


%% Analysis
t = linspace(0,(T-1)/Fs,T);
nu0 = Fs/4;
freqdisp = [16 8 4 2]; % Displayed frequencies in kHz
sdisp = log2(nu0./(1e3*freqdisp)); % coreesponding log-scales

figure;
for n=1:N
    Whyr = cwt_JEFAS(haty(n,:),scales,wav_typ,wav_param);
    subplot(N,1,n); imagesc(t,log2(scales),abs(Whyr));
    yticks(sdisp); yticklabels(freqdisp);
    xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray));
end
title('Estimated sources');

figure;
for n=1:N
    Whyr = cwt_JEFAS(z(n,:),scales,wav_typ,wav_param);
    subplot(N,1,n); imagesc(t,log2(scales),abs(Whyr));
    yticks(sdisp); yticklabels(freqdisp);
    xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray));
end
title('Observations')

% save('res_bss_nonstat','y','pileBoptim','dgammaML','Sx');
%% Performances
[SDR,SIR,SAR,perm] = bss_eval_sources(haty,y);

% Amari index
for k=1:Kmat
   indJEFAS(k) = amari(pileBoptim(:,:,k),pileA(:,:,vectau(k)));
   indSOBI(k) = amari(Bsobi,pileA(:,:,vectau(k)));
   indTFQ(k) = amari(BTFQ,pileA(:,:,vectau(k)));
   k2 = length(vectauns(vectauns<=vectau(k)));
   indPSOBI(k) = amari(pile_Bsobins(:,:,k2),pileA(:,:,vectau(k)));
end

fprintf('Critère | SOBI    | p-SOBI  | TFQ  | JEFAS-BSS \n')
fprintf('SIR     |  %.2f  |   %.2f  |  %.2f  |  %.2f\n', mean(SIRsobi),mean(SIRsobins),mean(SIRmoreau),mean(SIR))
fprintf('SDR     |  %.2f  |  %.2f  |  %.2f  |  %.2f\n', mean(SDRsobi),mean(SDRsobins),mean(SDRmoreau),mean(SDR))
fprintf('SAR     |  %.2f  |  %.2f  | %.2f  |  %.2f\n', mean(SARsobi),mean(SARsobins),mean(SARmoreau),mean(SAR))
fprintf('Amari   | %.2f  | %.2f  | %.2f  | %.2f\n', mean(indSOBI),mean(indPSOBI),mean(indTFQ),mean(indJEFAS))

figure;subplot(2,1,2);
plot(t(vectau),indJEFAS,t(vectau),indSOBI,'k--',t(vectau),indTFQ,'g-.',t(vectau),indPSOBI,'r:','linewidth',2); grid on; axis tight;
xlabel('Time (s)'); ylabel('Amari index (dB)'); grid on; legend('JEFAS-BSS','SOBI','TFQ','p-SOBI');
%% Convergence
iter = 1:length(errSAR);

subplot(2,1,1);
plot(iter,estSIR,'linewidth',2);
xlabel('Iteration'); ylabel('Convergence criterion (dB)'); grid on; axis tight;
