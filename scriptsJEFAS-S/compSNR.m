%% JEFAS-S initialized with the results of JEFAS

clear all; close all; clc;
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-S'));
addpath('../signals');
addpath('../cwt');

load('toysigJEFASS')
theta = log2(dgamma);
T = length(y0);
%% JEFAS-S

Nbscales = 100;
scales = 2.^linspace(0,3,Nbscales);
scales = scales(1:7:end); % pour se comparer a JEFAS
wav_param = 300;
wav_paramWP = 10; % used by JEFAS and JEFASS
wav_typ = 'sharp';

% parameters specified for initialization by JEFAS algorithm
ratio = 0.05; % activity threshold
dgamma0 = ones(1,T);
a0 = ones(1,T);
stop_crit = 5e-3; % relative update threshold
paramWAV = {wav_typ,wav_param,wav_paramWP};
paramWP = {scales};
paramAM = {'no AM'};
scalesS = 2.^(linspace(-1,3.5,150)); % for spectrum estimation
paramS = {scalesS};

itD = 10; % nb d'iterations de la descente de gradient
stopD = 1e-4; % tolerance de la descente de gradient
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxIterations',itD,'StepTolerance',stopD,'Display','off'); % fmincon

Dt = 50; % valeur du sous echantillonage pour l'estimation de thetaEM
TT = 512; %1024;%T/16;%256 % taille de la decoupe pour W
Delta = 0.75*TT; % recouvrement
alpha = 15; % largeur pour l'estimation du spectre par Welch

Nit = 5; % number of iterations
thres = 0.2; % ?

wav_paramR = 200;
scalesR = 2.^linspace(0,3,Nbscales);
[M_psi,M_tmpdpsi] = bas_calc_dcov(scalesR,wav_typ,wav_paramR,T);
MatPsi = ifft(M_psi.',[],1);

P = 30;
sigmaH = logspace(-3,-0.5,P);

for p = 1:P
    
    sigmay = sigmaH(p) ;
    y = y0 + sigmay*randn(T,1) ;
    
    % Apply JEFAS Initialization
    [~, dgammaML, ~, ~] = estim_altern(y,Dt,ratio,dgamma0,a0,paramWAV,paramWP,paramAM,paramS,stop_crit,10,0);
    thetaINIT = log2(dgammaML); % Initialization from JEFAS results
    
    % Apply JEFAS-S
    tic;
    [dgammaEST, SxEST, W, nll] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_paramWP,Dt,TT,Delta,alpha,Nit,thres,itD,stopD,theta);
    toc;

%     dgammaH{p} = dgammaEST ;
%     SxESTH{p} = SxEST ;
%     WH{p} = W ;

    uEST = abs(dgammaEST(:)-dgamma).^2;
    uML = abs(dgammaML(:)-dgamma).^2;

    fprintf('Noise variance %.2e\n\n', sigmaH(p)^2);
    fprintf('EQM:\n JEFAS: %.3e\n JEFAS-S: %.3e\n\n',mean(uML(1:8151)),mean(uEST(1:8151)));  % prevent edge effect
    
    [Wfinest, MMSigmay] = transform_adap(y,sigmay,TT,Delta,M_psi,SxEST,log2(dgammaEST),MatPsi); % adpated representation
    Cy = buildSigmay(MMSigmay,T,TT,Delta);
    yr = synthesis(Wfinest, MatPsi); % reconstruction

    snrINth = 10*log10(1 + var(y0)/sigmay^2);
    snrIN = snr(y,y-y0);
    snrOUT = snr(yr(:),yr(:)-y0);
    fprintf('Measured input SNR %.2f\nMeasured output SNR %.2f\n',snrIN,snrOUT)
    
    biasTH = theo_bias(y0,TT,Delta,MMSigmay,sigmay);
    biasXP = yr(:)-y0;

    % Theoretical signal cov matrix
    [~, MMSigmayTH] = transform_adap(y,sigmay,TT,Delta,M_psi,Sx,log2(dgamma),MatPsi); % adpated representation
    CyTH = buildSigmay(MMSigmayTH,T,TT,Delta);
    CyiTH = inv(CyTH);
    % Theoretical SNR
    varty0 = trace( (eye(T)-sigmay^2*CyiTH)^2 * (sigmay^2*eye(T) + y0*y0.') );
    vartnoise = sigmay^2 * ( trace((eye(T)-sigmay^2*CyiTH)^2) + trace(sigmay^2*CyiTH^2*y0*y0.') );
    snrOUTth = 10*log10( varty0 / vartnoise ); % Theoretical output SNR
    fprintf('Theoretical Output SNR %.2f\n\n',snrOUTth)

    fprintf('------------------------------------------------\n\n');
    
    snrinH(p) = snrIN ;
    snroutH(p) = snrOUT ;
    snroutthH(p) = snrOUTth ;
end

save('results/resultsSNR','snrinH','snroutH','snroutthH') ;
