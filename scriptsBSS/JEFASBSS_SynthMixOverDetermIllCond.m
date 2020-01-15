clear all; close all; clc;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-BSS'));

load('../signals/synthetic_NonstatMixtureCrossingOverDetermIllCond.mat');
Fs = 44100;

[Nz,T] = size(z);


%% Source number estimation
thres = 0.5;
winlength = round(logspace(-3.5, -0.3, 10)*T);
k = 1;
for Delta=winlength
    [~,vecSV] = estimSourcesNb(z,Delta,thres);
    sv = mean(vecSV,2);
    rap(k) = sv(Nz)/sv(Nz-1);
    k = k+1;
end

semilogx(winlength/Fs,rap,'-bo','linewidth',2); axis tight; grid on;
xlabel('Window Length (s)'); ylabel('Singular values ratio'); set(gca,'FontSize',18)

Ny = estimSourcesNb(z,round(0.01*T),thres);
%% JEFAS-BSS estimation

dgamma0 = ones(Ny,T);

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
paramS = {scalesS};

stop_crit = 5e-3;
stopSIR = 75; % stopping criterion

combinE = nchoosek(1:Nz,Ny);
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

% Reorder estimated sources and compute the global separation matrix
heapBtot = zeros(Ny,Nz,Kmat);
for indBSS = 1:nbBSS
    [hy{indBSS}, st{indBSS}] = reordersig(y,hy{indBSS});
    
    combin = combinE(indBSS,:);
    heapBoptim = hB{indBSS};
    for k=1:Kmat
       heapBtot(:,combin,k) = heapBtot(:,combin,k) + heapBoptim(st{indBSS},:,k);
    end
end
heapBtot = heapBtot/nchoosek(Nz-1,Ny-1);
hytot = nonstatunmixing(z,heapBtot,vectau); % final sepration source

dgammaprec = dg{1};
aprec = ones(Ny,T);
for n=1:Ny
    ytot = hytot(n,:); % estimated source n
    [aML,dgammaML, Sxn] = estim_altern(ytot,Dt,0,dgammaprec(n,:),aprec(n,:),paramWAV,paramWP,paramAM,paramS,stop_crit,10,1);
    dgammatot(n,:) = dgammaML;
    atot(n,:) = aML;
    Sxtot(n,:) = Sxn;
end

%% Performance evaluation
close all;
t = linspace(0,(T-1)/Fs,T);

 for n=1:T
    A = heapA(:,:,n);
    c1(n)=cond(A([1,2],:));
    c2(n)=cond(A([1,3],:));
    c3(n)=cond(A([2,3],:));
end
figure;
subplot(3,1,1); plot(t,c1); xlabel('Time (s)'); ylabel('Condition number'); title('Submatrix #1');
subplot(3,1,2); plot(t,c2); xlabel('Time (s)'); ylabel('Condition number'); title('Submatrix #2');
subplot(3,1,3); plot(t,c3); xlabel('Time (s)'); ylabel('Condition number'); title('Submatrix #3');


figure;
for indBSS = 1:nbBSS
    combin = combinE(indBSS,:);
    zcom = z(combin,:);
    
    heapBoptim = hB{indBSS};
    haty = hy{indBSS}; % unmixing
    
    [SDR{indBSS},SIR{indBSS}] = bss_eval_sources(haty,y); % JEFAS-BSS
    for k=1:Kmat
       indJEFAS(k) = amari(heapBoptim(:,:,k),heapA(combin,:,vectau(k)));
    end
    indam{indBSS} = indJEFAS;
    
    hold on; plot(t(vectau),indJEFAS,'linewidth',2,'DisplayName',['Combination #' num2str(indBSS)]); grid on; axis tight; 
    xlabel('Time (s)'); ylabel('Amari index (dB)'); grid on; set(gca,'fontsize',24);
end
legend show; xlabel('Time (s)'); ylabel('Amari index (dB)');


[SDRtot,SIRtot] = bss_eval_sources(hytot,y); % JEFAS-BSS

for k=1:Kmat
   indJEFAStot(k) = amari(heapBtot(:,:,k),heapA(:,:,vectau(k)));
end

fprintf('Index   | JEFAS-BSS | s. JEFAS-BSS #1 | s. JEFAS-BSS #2 | s. JEFAS-BSS #3 \n')
fprintf('SIR     |   %.2f   |      %.2f      |      %.2f      |      %.2f\n', mean(SIRtot), mean(SIR{1}), mean(SIR{2}), mean(SIR{3}) )
fprintf('SDR     |   %.2f   |      %.2f      |      %.2f      |       %.2f\n', mean(SDRtot), mean(SDR{1}), mean(SDR{2}), mean(SDR{3}) )
fprintf('Amari   |  %.2f   |     %.2f      |     %.2f      |     %.2f\n', mean(indJEFAStot), mean(indam{1}), mean(indam{1}), mean(indam{3}) )


%% OLD
% nu0 = Fs/4;
% freqdisp = [8 4 2 1 0.5]; % Displayed frequencies in kHz
% sdisp = log2(nu0./(1e3*freqdisp)); % corresponding log-scales
% 
% figure; title('Original Sources and observations');
% for n=1:Nz
%     if n <= Ny
%         Wy = cwt_JEFAS(y(n,:),scales,wav_typ,wav_param);
%         subplot(2,Nz,n); imagesc(t,log2(scales),abs(Wy)); yticks(sdisp); yticklabels(freqdisp);
%         ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
%     end
%     Wz = cwt_JEFAS(z(n,:),scales,wav_typ,wav_param);
%     subplot(2,nbBSS,Nz+n); imagesc(t,log2(scales),abs(Wz)); yticks(sdisp); yticklabels(freqdisp);
%     xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
% end
% 
% figure; title('JEFAS-BSS estimates of the sources');
% for n=1:Ny
%     
%     Why = cwt_JEFAS(hytot(n,:),scales,wav_typ,wav_param);
%     subplot(3,Ny,n); imagesc(t,log2(scales),abs(Why));
%     yticks(sdisp); yticklabels(freqdisp);
%     xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
%     
%     subplot(3,Ny,Ny+n); plot(t,dgamma(n,:),'b--',t,dgammatot(n,:),'r','linewidth',2);
%     xlabel('Time (s)'); ylabel('\gamma''(t)'); grid on; legend('Ground truth function','Estimated function'); set(gca,'fontsize',18);
%     
%     hatx = statAMWP(hytot(n,:),atot(n,:),dgammatot(n,:));
%     alpha = 15;
%     Nff = 50000;
%     Sxw = estim_spec(hatx,Nff,alpha);
%     Sxw = Sxw*sum(Sx(n,:)  + 25e-4)/sum(Sxw); % normalization
%     freq = linspace(0,Fs,Nff);
%     freq2 = linspace(0,(T-1)*Fs/T,Fs);
%     subplot(3,Ny,2*Ny+n); plot(freq2/1e3,log10(Sx(n,:)  + 25e-4),'b--',freq/1e3,log10(Sxw),'r','linewidth',2); 
%     xlabel('Frequency (kHz)'); ylabel('S_x'); axis tight; grid on;
%     xlim([0 13]); ylim([-5.1 1]); yticks([-5 -4 -3 -2 -1 0]); yticklabels({'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '10^{0}'}); 
% end