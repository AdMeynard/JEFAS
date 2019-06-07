clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-BSS'));

load('../signals/synthetic_NonstatMixtures.mat'); 

Fs = 44100;

K = 20; % Number of synthetic signals
for l = 1:K
    close all;
    fprintf('SYNTHETIC SIGNAL BSS # %i  \n\n', l);
    
    y = yK(:,:,l);
    z = zK(:,:,l);
    heapA = heapAK(:,:,:,l);
    heapB = heapBK (:,:,:,l);
    Sx = SxK (:,:,l);
    dgamma = dgammaK (:,:,l); 
    
    
    [N,T] = size(z);

    Kmat = 200; % number of instants where we estimate the unmixing matrix
    vectau = floor(linspace(1,T-1,Kmat));

    %% SOBI estimation
    p = 1000; % Number of jointly diagonalized matrices
    [Asobi,ysobi] = sobi(z,N,p);
    Bsobi = inv(Asobi);

    %% p-SOBI estimation
    Dtp = 4000;
    vectauns = 1:Dtp:T;
    [heap_Apsobi, heap_Bpsobi] = psobi(z,vectauns);
    haty0 = nonstatunmixing(z,heap_Bpsobi,vectauns);

    %% QTF-BSS estimation
    eps3 = 0.1; % imaginary part threshold
    eps4 = 100; % real part threshold

    pp = 10; % subsampling
    nn = N; % number of classes
    AQTF = BSS_QTF(z, pp, eps3, eps4, nn); % QTF BSS
    BQTF = inv(AQTF);

    hatyqtf = BQTF*z;

    %% JEFAS-BSS estimation
    dgamma0 = ones(N,T);

    Dt = 200;

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

    paramBSS = {scales,vectau,L_bss};
    paramWAV = {wav_typ,wav_param,wav_paramWP};
    paramWP = {scalesWP,0.05};
    paramAM = {'no AM'};
%     paramAM = {'AM',scales,rAM};
    paramS = {scalesS};

    stop_crit = 5e-3;
    stopSIR = 75; % stopping criterion

    init_meth = 'sobi';
    [heapB_init, vectau_init] = JEFASBSSinit(z, init_meth, Dtp);

    tic;
    [heapBoptim,dgammaML,aML,SxML,convergeSIR] = altern_bss_deform_nonstat(z,heapB_init,vectau_init,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,stopSIR);
    tK(l) = toc;

    haty = nonstatunmixing(z,heapBoptim,vectau); % unmixing

    %% Performances
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
    
    SIRsobiK(l) = mean(SIRsobi);
    SIRpsobiK(l) = mean(SIRpsobi);
    SIRqtfK(l) = mean(SIRqtf);
    SIRK(l) = mean(SIR);
    SDRsobiK(l) = mean(SDRsobi);
    SDRpsobiK(l) = mean(SDRpsobi);
    SDRqtfK(l) = mean(SDRqtf);
    SDRK(l) = mean(SDR);
    indSOBIK(l) = mean(indSOBI);
    indPSOBIK(l) = mean(indPSOBI);
    indQTFK(l) = mean(indQTF);
    indJEFASK(l) = mean(indJEFAS);


end

fprintf('Algorithm ||      SDR      ||      SIR     ||  Amari index \n')
fprintf('          ||  Mean  |  SD  ||  Mean |  SD  ||  Mean |  SD   \n')
fprintf('SOBI      ||  %.2f  | %.2f || %.2f | %.2f || %.2f | %.2f\n', mean(SDRsobiK),std(SDRsobiK),mean(SIRsobiK),std(SIRsobiK),mean(indSOBIK),std(indSOBIK))
fprintf('p-SOBI    ||  %.2f | %.2f ||  %.2f | %.2f || %.2f | %.2f\n', mean(SDRpsobiK),std(SDRpsobiK),mean(SIRpsobiK),std(SIRpsobiK),mean(indPSOBIK),std(indPSOBIK))
fprintf('QTF-BSS   ||  %.2f | %.2f ||  %.2f | %.2f || %.2f | %.2f\n', mean(SDRqtfK),std(SDRqtfK),mean(SIRqtfK),std(SIRqtfK),mean(indQTFK),std(indQTFK))
fprintf('JEFAS-BSS ||  %.2f | %.2f || %.2f | %.2f ||%.2f | %.2f\n', mean(SDRK),std(SDRK),mean(SIRK),std(SIRK),mean(indJEFASK),std(indJEFASK))


