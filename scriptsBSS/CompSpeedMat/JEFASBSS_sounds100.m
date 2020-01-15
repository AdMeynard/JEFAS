clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-BSS'));

load('../signals/soundMixture100.mat');
[N,T] = size(z);

%% SOBI estimation
p = 1000; % Number of jointly diagonalized matrices
[Asobi,ysobi] = sobi(z,N,p);
Bsobi = inv(Asobi);

%% p-SOBI estimation
Dtp = 4000;
vectauns = 1:Dtp:T;
[heap_Apsobi, heap_Bpsobi] = psobi(z,N,vectauns);
haty0 = nonstatunmixing(z,heap_Bpsobi,vectauns);

%% QTF-BSS estimation
eps3 = 0.1; % imaginary part threshold
eps4 = 100; % real part threshold

pp = 100; % subsampling
nn = N; % number of classes
AQTF = BSS_QTF(z, pp, eps3, eps4, nn); % QTF BSS
BQTF = inv(AQTF);

hatyqtf = BQTF*z;

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
t = linspace(0,(T-1)/Fs,T);
nu0 = Fs/4;
freqdisp = [16 8 4 2 1 0.5 0.25]; % Displayed frequencies in kHz
sdisp = log2(nu0./(1e3*freqdisp)); % corresponding log-scales

[haty, stablematch] = reordersig(y,haty);

%% Performances
[SDR0,SIR0] = bss_eval_sources(z,y); % no unmixing
[SDRsobi,SIRsobi] = bss_eval_sources(ysobi,y); % SOBI
[SDRpsobi,SIRpsobi] = bss_eval_sources(haty0,y); %p-SOBI
[SDRqtf,SIRqtf] = bss_eval_sources(hatyqtf,y); % QTF-BSS
[SDR,SIR] = bss_eval_sources(haty,y); % JEFAS-BSS

% Amari index
for k=1:Kmat
   indJEFAS(k) = amari(heapBoptim(:,:,k),heapA(:,:,vectau(k)));
   indSOBI(k) = amari(Bsobi,heapA(:,:,vectau(k)));
   indQTF(k) = amari(BQTF,heapA(:,:,vectau(k)));
   k2 = length(vectauns(vectauns<=vectau(k)));
   indPSOBI(k) = amari(heap_Bpsobi(:,:,k2),heapA(:,:,vectau(k)));
end

%% Save
save('../results/results100','indJEFAS','indSOBI','indQTF','indPSOBI','SDR0','SIR0','SDRsobi','SIRsobi','SDRpsobi','SIRpsobi','SDRqtf','SIRqtf','SDR','SIR');