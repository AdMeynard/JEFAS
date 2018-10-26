clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath('JEFAS-BSS');
addpath(genpath('othermethods_BSS')); % quadratic TF BSS method

% load('../signals/sig_compBSS.mat'); Fs = 44100; heapA = repmat(A,1,1,Fs);
load('../signals/sig_compBSS_nonstatK.mat'); Fs = 44100;

K = 5;
for l = 1:K
    y = yK(:,:,l);
    z = zK(:,:,l);
    heapA = heapAK(:,:,:,l);
    heapB = heapBK (:,:,:,l);
    Sx = SxK (:,:,l);
    dgamma = dgammaK (:,:,l); 
    
    
    [N,T] = size(z);

    Kmat = 200; % nombre d'instants ou l'on calcule Boptim
    vectau = floor(linspace(1,T-1,Kmat));

    %% 1ere methode: via SOBI
    p = 1000; % Nombre de matrices à diagonaliser conjointement
    [Asobi,ysobi] = sobi(z,N,p);
    Bsobi = inv(Asobi);
    Bsobi = Bsobi./sqrt(sum(abs(Bsobi).^2,2));

    [SDRsobi,SIRsobi,SARsobi,permsobi] = bss_eval_sources(ysobi,y); % SOBI

    %% 2eme methode: via p-SOBI
    Dtp = 4000;
    vectauns = 1:Dtp:T;
    [heap_Apsobi, heap_Bpsobi] = psobi(z,vectauns);
    haty0 = nonstatunmixing(z,heap_Bpsobi,vectauns);

    [SDRpsobi,SIRpsobi,SARpsobi,permpsobi] = bss_eval_sources(haty0,y);

    %% 3eme methode: via QTF
    eps3 = 0.1;
    eps4 = 100;

    pp = 10; % subsampling
    nn = 2; % number of classes
    AQTF = BSS_QTF(z, pp, eps3, eps4,nn); % QTF BSS
    BQTF = inv(AQTF);

    hatyqtf = BQTF*z;
    [SDRqtf,SIRqtf,SARqtf,permqtf] = bss_eval_sources(hatyqtf,y);

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
    cvmax = 15;
    stopSIR = 75; % critère d'arret sur l'amélioration

    init_meth = 'psobi';
    [heapB_init, vectau_init] = JEFASBSSinit(z, init_meth, Dtp);

    tic;
    [heapBoptim,dgammaML,aML,SxML,errSIR,errSAR,errSDR,estSIR] = altern_bss_deform_nonstat(y,z,heapB_init,vectau_init,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,cvmax,stopSIR,options);
    toc;

    haty = nonstatunmixing(z,heapBoptim,vectau); % unmixing


    % save('res_bss_nonstat','y','heapBoptim','dgammaML','Sx');
    %% Performances
    [SDR,SIR,SAR,perm] = bss_eval_sources(haty,y);

    % Amari index
    for k=1:Kmat
       indJEFAS(k) = amari(heapBoptim(:,:,k),heapA(:,:,vectau(k)));
       indSOBI(k) = amari(Bsobi,heapA(:,:,vectau(k)));
       indQTF(k) = amari(BQTF,heapA(:,:,vectau(k)));
       k2 = length(vectauns(vectauns<=vectau(k)));
       indPSOBI(k) = amari(heap_Bpsobi(:,:,k2),heapA(:,:,vectau(k)));
    end

    fprintf('Critère | SOBI    | p-SOBI  | QTF  | JEFAS-BSS \n')
    fprintf('SIR     |  %.2f  |   %.2f  |  %.2f  |  %.2f\n', mean(SIRsobi),mean(SIRpsobi),mean(SIRqtf),mean(SIR))
    fprintf('SDR     |  %.2f  |  %.2f  |  %.2f  |  %.2f\n', mean(SDRsobi),mean(SDRpsobi),mean(SDRqtf),mean(SDR))
    fprintf('SAR     |  %.2f  |  %.2f  | %.2f  |  %.2f\n', mean(SARsobi),mean(SARpsobi),mean(SARqtf),mean(SAR))
    fprintf('Amari   | %.2f  | %.2f  | %.2f  | %.2f\n', mean(indSOBI),mean(indPSOBI),mean(indQTF),mean(indJEFAS))
    
    SIRsobiK(l) = mean(SIRsobi);
    SIRpsobiK(l) = mean(SIRpsobi);
    SIRqtfK(l) = mean(SIRqtf);
    SIRK(l) = mean(SIR);
    indSOBIK(l) = mean(indSOBI);
    indPSOBIK(l) = mean(indPSOBI);
    indQTFK(l) = mean(indQTF);
    indJEFASK(l) = mean(indJEFAS);


end
