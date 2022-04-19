%% JEFAS-S on the narrowband signal

clear all; close all; clc;
addpath(genpath('../JEFASalgo'));
addpath(genpath('../JEFAS-S'));
addpath('../signals');
addpath('../cwt');
addpath('results');

cd ../JEFAS-S/minFunc_2012       % Change to the unzipped directory
mexAll;               % Compile mex files (not necessary on all systems)

cd ../../scriptsJEFAS-S

load('TwoSineWaves_fast')
theta = log2(dgamma);
T = length(y);
%% JEFAS-S

Nbscales = 125 ;
scales = 2.^linspace(-0.67,2.33,Nbscales);
subrate = 7 ; % subsampling step for the scales to ensure the covariance invertibility
scales = scales(1:subrate:end); % for comparison with JEFAS
wav_param = 20;
wav_typ = 'sharp';

prior = 'wavelet' ; %prior for the covariance matrix of the time scale representation
priorList = {prior} ;

itD = 10; % number of iterations for gradient descent
stopD = 1e-4; % gradient descent: tolerance
options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'MaxIterations',itD,'StepTolerance',stopD,'Display','off'); % fmincon

Dt = 1; % subsampling value for thetaEM estimation
TT = 512; %1024;%T/16;%256 %  slicing size for W
Delta = 0.75*TT; % overlap
alpha = 7 ; % width for Welch estimation

TS = 1024 ;

Nit = 300 ; % number of iterations

initMeth = 'none';
load('JEFAS_2sinesfast.mat')
switch initMeth
    
    case 'true'
        thetaINIT = log2(dgamma);
    
    case 'none'
        thetaINIT = zeros(1,T);
        
    case 'JEFAS'
        thetaINIT = log2(dgammaML); % Initialization from JEFAS results

end

thres = 0.5;
tic;
[dgammaEST, SxEST, W, nll] = EMwarping(y,sigmay,thetaINIT,scales,wav_typ,wav_param,priorList,Dt,TT,Delta,TS,alpha,Nit,thres,itD,stopD,theta);
toc;

save('results/JEFASS_2sinesfast','W', 'dgammaEST','SxEST','y','Fs')

%% Results
t = linspace(0,(T-1)/Fs,T);

figure;
plot(nll,'linewidth',2); grid on;
xlabel('Iteration'); ylabel('negative log-likelihood of the signal');

figure;
plot(t,dgamma,'k--',t,dgammaML,'b',t,dgammaEST,'r','linewidth',2);
xlabel('Time (s)'); ylabel("\gamma'");
legend('Ground truth function','JEFAS estimate','JEFAS-S estimate');
axis tight; grid on; ylim([0.5 1.5]);
set(gca,'fontsize',18);

uEST = abs(dgammaEST(:)-dgamma(:)).^2;
uML = abs(dgammaML(:)-dgamma(:)).^2;
fprintf('MSE:\n JEFAS: %.3e\n JEFAS-S: %.3e\n',mean(uML),mean(uEST)) ;
