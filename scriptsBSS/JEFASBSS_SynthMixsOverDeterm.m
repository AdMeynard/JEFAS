clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-BSS'));

load('../signals/synthetic_NonstatMixturesOverDet.mat'); 

Fs = 44100;

K = 20; % Number of synthetic signals
thres = 0.5; %Threshold for source number estimation
for l = 1:K
    close all;
    fprintf('SYNTHETIC SIGNAL BSS # %i  \n\n', l);
    
    y = yK(:,:,l);
    z = zK(:,:,l);
    heapA = heapAK(:,:,:,l);
    Sx = SxK (:,:,l);
    dgamma = dgammaK (:,:,l); 
    
    
    [Nz,T] = size(z);
    Ny = estimSourcesNb(z,round(0.01*T),thres);

    Kmat = 200; % number of instants where we estimate the unmixing matrix
    vectau = floor(linspace(1,T-1,Kmat));

    %% SOBI estimation
    p = 1000; % Number of jointly diagonalized matrices
    

    %% p-SOBI estimation
    Dtp = 4000;
    vectauns = 1:Dtp:T;

    %% QTF-BSS estimation
    eps3 = 0.1; % imaginary part threshold
    eps4 = 100; % real part threshold

    pp = 10; % subsampling
    nn = Ny; % number of classes
    
    combinE = nchoosek(1:Nz,Ny);
    nbBSS = size(combinE,1);
    
    BQTF = zeros(Ny,Nz);
    Bsobi = zeros(Ny,Nz);
    heap_Bpsobi = zeros(Ny,Nz,length(vectauns));
    for indBSS = 1:nbBSS
        combin = combinE(indBSS,:);
        zcom = z(combin,:);
        
        [Asobicom,ysobi] = sobi(zcom,p);
        Bsobicom = pinv(Asobicom);
        hatysobicom = Bsobicom*zcom;
        [~, st] = reordersig(y,hatysobicom);
        Bsobi(:,combin) = Bsobi(:,combin) + Bsobicom(st,:);
        
        [heap_Apsobicom, heap_Bpsobicom] = psobi(zcom,Ny,vectauns);
        haty0 = nonstatunmixing(zcom,heap_Bpsobicom,vectauns);
        [~, st] = reordersig(y,haty0);
        for k=1:length(vectauns)
            heap_Bpsobi(:,combin,k) = heap_Bpsobi(:,combin,k) + heap_Bpsobicom(st,:,k);
        end
        
        AQTFcom = BSS_QTF(zcom, pp, eps3, eps4, nn); % QTF BSS
        BQTFcom = pinv(AQTFcom);
        hatycom = BQTFcom*zcom;
        [~, st] = reordersig(y,hatycom);
        BQTF(:,combin) = BQTF(:,combin) + BQTFcom(st,:);
    end
    
    Bsobi = Bsobi/nchoosek(Nz-1,Ny-1);
    hatysoby = Bsobi*z;
    
    heap_Bpsobi = heap_Bpsobi/nchoosek(Nz-1,Ny-1);
    haty0 = nonstatunmixing(z,heap_Bpsobi,vectauns);
    
    BQTF = BQTF/nchoosek(Nz-1,Ny-1);
    hatyqtf = BQTF*z;

    %% JEFAS-BSS estimation
    dgamma0 = ones(Ny,T);

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
    paramS = {scalesS};

    stop_crit = 5e-3;
    stopSIR = 75; % stopping criterion

    init_meth = 'sobi';
    heapBtot = zeros(Ny,Nz,Kmat);
    for indBSS = 1:nbBSS
        combin = combinE(indBSS,:);
        zcom = z(combin,:);
        [heapB_init, vectau_init] = JEFASBSSinit(zcom, init_meth, Dtp);

        tic;
        [heapBoptim,dgammaML,aML,SxML,convergeSIR] = altern_bss_deform_nonstat(zcom,heapB_init,vectau_init,dgamma0,Dt,paramBSS,paramWAV,paramWP,paramAM,paramS,stop_crit,stopSIR);
        toc;
        hB{l,indBSS} = heapBoptim;
        dg{l,indBSS} = dgammaML;
        SxE{l,indBSS} = SxML;
        cS{l,indBSS} = convergeSIR;

        haty = nonstatunmixing(zcom,heapBoptim,vectau); % unmixing
        [haty, st] = reordersig(y,haty);
        for k=1:Kmat
           heapBtot(:,combin,k) = heapBtot(:,combin,k) + heapBoptim(st,:,k);
        end
        
        %performances
        [SDRs{l,indBSS},SIRs{l,indBSS}] = bss_eval_sources(haty,y); % JEFAS-BSS
        for k=1:Kmat
           indJEFASs(k) = amari(heapBoptim(:,:,k),heapA(combin,:,vectau(k)));
        end
        indamJEFASs{l,indBSS} = indJEFASs;
    end
    
    heapBtot = heapBtot/nchoosek(Nz-1,Ny-1);
    hatytot = nonstatunmixing(z,heapBtot,vectau); % final sepration source
    
    hBtot{l} = heapBtot;
    hytot{l} = hatytot;
   % Reorder estimated sources and compute the global separation matrix
%     dgammaprec = dg{l,1};
%     aprec = ones(Ny,T);
%     for n=1:Ny
%         ytot = hytot(n,:); % estimated source n
%         [aML,dgammaML, Sxn] = estim_altern(ytot,Dt,0,dgammaprec(n,:),aprec(n,:),paramWAV,paramWP,paramAM,paramS,stop_crit,10,1);
%         dgammatot(n,:) = dgammaML;
%         atot(n,:) = aML;
%         Sxtot(n,:) = Sxn;
%     end

    %% Performances
    [SDRsobi,SIRsobi] = bss_eval_sources(ysobi,y); % SOBI
    [SDRpsobi,SIRpsobi] = bss_eval_sources(haty0,y); %p-SOBI
    [SDRqtf,SIRqtf] = bss_eval_sources(hatyqtf,y); % QTF-BSS
    [SDR,SIR] = bss_eval_sources(hatytot,y); % JEFAS-BSS
    for indBSS = 1:nbBSS
        SDRss(indBSS) = mean(SDRs{l,indBSS}); % JEFAS-BSS single
        SIRss(indBSS) = mean(SIRs{l,indBSS}); % JEFAS-BSS single
    end
    SDRss = mean(SDRss);
    SIRss = mean(SIRss);

    % Amari index
    for k=1:Kmat
       indJEFASBSS(k) = amari(heapBtot(:,:,k),heapA(:,:,vectau(k)));
       indSOBI(k) = amari(Bsobi,heapA(:,:,vectau(k)));
       indQTF(k) = amari(BQTF,heapA(:,:,vectau(k)));
       k2 = length(vectauns(vectauns<=vectau(k)));
       indPSOBI(k) = amari(heap_Bpsobi(:,:,k2),heapA(:,:,vectau(k)));
    end
    indJEFASBSSs = zeros(1,Kmat);
    for indBSS = 1:nbBSS
        indJEFASBSSs = indJEFASBSSs + indamJEFASs{l,indBSS};
    end
    indJEFASBSSs = indJEFASBSSs/nbBSS;
    
    SIRsobiK(l) = mean(SIRsobi);
    SIRpsobiK(l) = mean(SIRpsobi);
    SIRqtfK(l) = mean(SIRqtf);
    SIRK(l) = mean(SIR);
    SIRsK(l) = mean(SIRss);
    SDRsobiK(l) = mean(SDRsobi);
    SDRpsobiK(l) = mean(SDRpsobi);
    SDRqtfK(l) = mean(SDRqtf);
    SDRK(l) = mean(SDR);
    SDRsK(l) = mean(SDRss);
    indSOBIK(l,:) = indSOBI;
    indPSOBIK(l,:) = indPSOBI;
    indQTFK(l,:) = indQTF;
    indJEFASK(l,:) = indJEFASBSS;
    indJEFASBSSsK(l,:) = indJEFASBSSs;
end

fprintf('Algorithm        ||      SDR      ||      SIR     ||  Amari index \n')
fprintf('                 ||  Mean  |  SD  ||  Mean |  SD  ||  Mean |  SD   \n')
fprintf('SOBI             ||   %.2f | %.2f || %.2f | %.2f || %.2f | %.2f\n', mean(SDRsobiK),std(SDRsobiK),mean(SIRsobiK),std(SIRsobiK),mean(indSOBIK(:)),std(indSOBIK(:)))
fprintf('p-SOBI           ||  %.2f | %.2f ||  %.2f | %.2f || %.2f | %.2f\n', mean(SDRpsobiK),std(SDRpsobiK),mean(SIRpsobiK),std(SIRpsobiK),mean(indPSOBIK(:)),std(indPSOBIK(:)))
fprintf('QTF-BSS          ||  %.2f | %.2f || %.2f | %.2f || %.2f | %.2f\n', mean(SDRqtfK),std(SDRqtfK),mean(SIRqtfK),std(SIRqtfK),mean(indQTFK(:)),std(indQTFK(:)))
fprintf('Single JEFAS-BSS ||   %.2f | %.2f || %.2f | %.2f ||%.2f | %.2f\n', mean(SDRsK),std(SDRsK),mean(SIRsK),std(SIRsK),mean(indJEFASBSSsK(:)),std(indJEFASBSSsK(:)))
fprintf('JEFAS-BSS        ||  %.2f | %.2f || %.2f | %.2f ||%.2f | %.2f\n', mean(SDRK),std(SDRK),mean(SIRK),std(SIRK),mean(indJEFASK(:)),std(indJEFASK(:)))



indSOBIm = mean(indSOBIK);
indPSOBIm = mean(indPSOBIK);
indQTFm = mean(indQTFK);
indJEFASm = mean(indJEFASK);
indJEFASBSSsm = mean(indJEFASBSSsK);

t = linspace(0,(T-1)/Fs,T);
figure;
plot(t(vectau),indSOBIm,'r:',t(vectau),indPSOBIm,'m',t(vectau),indQTFm,'k-.',t(vectau),indJEFASBSSsm,'c:',t(vectau),indJEFASm,'b','linewidth',2); grid on; axis tight; 
xlabel('Time (s)'); ylabel('Amari index (dB)'); grid on;
legend('SOBI','p-SOBI','QTF-BSS','Single JEFAS-BSS','JEFAS-BSS'); set(gca,'fontsize',21);
