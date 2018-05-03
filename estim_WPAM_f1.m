clear all; close all; clc;

warning off;
addpath('cwt');
addpath('deform_estimation');

load('signals/f1doppler');
T = length(y);
%% Joint estimation

Dt = 100; % temporal subsampling for the deformation estimation
dgamma0 = ones(1,T); % gamma'(t) initialization
a0 = ones(1,T); % a(t) initialization

wav_typ = 0; % wavelet type (cf. cwt.m)
wav_paramWP = 20; % corresponding parameter for warping estimation
wav_param = 500; % corresponding parameter for spectrum and AM estimations

NbScales = 125;
scalesAM = 2.^(linspace(1,6,NbScales));
subrate = 3; % subsampling step for the scales to ensure the covariance invertibility
scalesWP = scalesAM(1:subrate:end);

stopWP = 2e-2; % minimal gap between two steps in the gradient
itWP = 6; % number of gradient iterations

Nf = 2500; % number of frequencies for spectrum estimation
NbScalesS = 110;
scalesS = 2.^(linspace(-1,7,NbScalesS)); % for spectrum estimation

Nit = 10; % maximal number of iterations in the joint estimation
stop_crit = 5e-3; % relative update threshold

paramWAV = {wav_typ,wav_param,wav_paramWP};
paramWP = {scalesWP,itWP,stopWP};
paramS = {scalesS,Nf};

% WP estimation only
paramAM = {'no AM'}; % model with time warping only
tic;
[aML, dgammaML, Sx, evol_crit] = estim_altern(y,Dt,dgamma0,a0,paramWAV,paramWP,paramAM,paramS,stop_crit,Nit);
toc;
