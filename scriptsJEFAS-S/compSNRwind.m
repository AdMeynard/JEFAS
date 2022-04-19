%% Synthesis SNR as a function of the input noise 

clear all; close all; clc;
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-S'));
addpath('../signals');
addpath('../cwt');

load('wind4JEFASS');
y0 = y0(:) ;
theta = log2(dgamma);
T = length(y0);

Sx = sum(y0.^2)*Sx/sum(Sx);
%% JEFAS parameters

Dt = 2 ; % temporal subsampling for the deformation estimation
ratio = 0.05; % activity threshold
dgamma0 = ones(1,T); % gamma'(t) initialization
a0 = ones(1,T); % a(t) initialization

wav_param = 100;
wav_paramWP = 20;
wav_typ = 'sharp';

NbScales = 125;
scalesAM = 2.^(linspace(-0.75,2.12,NbScales));
subrate = 7 ; % subsampling step for the scales to ensure the covariance invertibility
scalesWP = scalesAM(1:subrate:end);

r = 1e-5; % regularization parameter

NbScalesS = 110;
scalesS = 2.^(linspace(-1,7,NbScalesS)); % for spectrum estimation

stop_crit = 1e-3; % relative update threshold

Nit = 30 ;

paramWAV = {wav_typ,wav_param,wav_paramWP};
paramWP = {scalesWP};
paramAM = {'no AM'}; % model without amplitude modulation
paramS = {scalesS};

%% Representation parameters

prior = 'wavelet' ; %prior for the covariance matrix of the time scale representation
priorList = {prior} ;

TT = T ;
Delta = 0.75*TT; % overlap
alpha = 7 ; % width for Welch estimation

[M_psi,M_tmpdpsi] = bas_calc_dcov(scalesAM,wav_typ,wav_param,T);
MatPsi = ifft(M_psi.',[],1);

%% SNR evolution

P = 30 ;
sigmaH = logspace(-1.67,0,P);

for p = 1:P
    
    sigmay = sigmaH(p) ;
    y = y0 + sigmay*randn(T,1) ;

    snrINth = 10*log10(var(y0)/sigmay^2);
    snrIN = snr(y0,y-y0);

    % Theoretical representation
    [WfinestTH, MMSigmayTH] = transform_adap(y,sigmay,priorList,TT,Delta,M_psi,Sx,log2(dgamma),MatPsi); % adapted representation
    CyTH = buildSigmay(MMSigmayTH,T,TT,Delta);
    yrTH = synthesis(WfinestTH, MatPsi); % reconstruction

    snrOUTexact = snr(y0(:),yrTH(:)-y0);
    fprintf('Exact Output SNR %.2f\n',snrOUTexact)

    % Theoretical SNR
    CyiTH = inv(CyTH) ;
    vartnoise = (4/T) * (sigmay^2 * ( trace((eye(T)-sigmay^2*CyiTH)^2) + trace(sigmay^2*CyiTH^2*(y0*y0')) )) ;
    snrOUTth = 10*log10( var(y0) / vartnoise ); % Theoretical output SNR
    fprintf('Theoretical Output SNR %.2f\n\n',snrOUTth)

    fprintf('------------------------------------------------\n\n');
    
    snrinH(p) = snrIN ;
    snroutexactH(p) = snrOUTexact ;
    snroutthH(p) = snrOUTth ;
end

save('results/resultsSNRwind','snrinH','snroutexactH','snroutthH') ;

figure;
plot(snrinH,snroutexactH,snrinH,snrinH,'k--','linewidth',2);
axis tight;
xlabel('Input SNR (dB)'); ylabel('Output SNR (dB)') ; grid on ;
set(gca,'fontsize',20);
