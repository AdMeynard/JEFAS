%% Spectral estimation from a dolphin sound recording
% Copyright (C) 2018 Adrien MEYNARD
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Author: Adrien MEYNARD
% Email: adrien.meynard@univ-amu.fr
% Created: 2018-05-23

%% Load signal
clear all; close all; clc;

warning off;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));

load('../signals/dolphin');
T = length(y);
%% Joint estimation

Dt = 100; % temporal subsampling for the deformation estimation
ratio = 0.05; % activity threshold
dgamma0 = ones(1,T); % gamma'(t) initialization
a0 = ones(1,T); % a(t) initialization

wav_typ = 'sharp'; % wavelet type (cf. cwt_JEFAS.m)
wav_paramWP = 20; % corresponding parameter for warping estimation
wav_param = 500; % corresponding parameter for spectrum and AM estimations

NbScales = 125;
scalesAM = 2.^(linspace(1,7,NbScales));
subrate = 3; % subsampling step for the scales to ensure the covariance invertibility
scalesWP = scalesAM(1:subrate:end);

r = 1e-5; % regularization parameter

NbScalesS = 110;
scalesS = 2.^(linspace(-1,7,NbScalesS)); % for spectrum estimation

stop_crit = 5e-3; % relative update threshold

paramWAV = {wav_typ,wav_param,wav_paramWP};
paramWP = {scalesWP};
paramS = {scalesS};

paramAM = {'AM',scalesAM,r};
tic;
[aML, dgammaML, Sx, evol_crit] = estim_altern(y,Dt,ratio,dgamma0,a0,paramWAV,paramWP,paramAM,paramS,stop_crit);
toc;


%% Analysis
t = 0:(1/Fs):((T-1)/Fs);

% Deformations:
figure;
subplot(1,2,1);plot(t,dgammaML,'linewidth',2); 
xlabel('Time (s)'); ylabel('Estimated \gamma''(t)'); axis tight; grid on;
subplot(1,2,2);plot(t,aML,'linewidth',2); 
xlabel('Time (s)'); ylabel('Estimated a^2(t)'); axis tight; grid on; ylim([0 2]);

% Spectrum:
z = statAMWP(y,aML,dgammaML);

alpha = 15;
Nff = 50000;
Sxw = estim_spec(z,Nff,alpha);
Syw = estim_spec(y,Nff,alpha);
freq = linspace(0,Fs,Nff);

figure;
subplot('Position', [0.06 0.55, 0.93, 0.45]);
plot(freq,log10(Syw),'linewidth',2); xlim([0 15000]); ylim([-6 1.5]); yticks([-6 -4 -2 0]); yticklabels({'10^{-6}' '10^{-4}' '10^{-2}' '10^{0}'});
ylabel('Estimated spectrum'); xticklabels([]); grid on;
set(gca,'fontsize',18);
subplot('Position', [0.06 0.08, 0.93, 0.45]);
plot(freq,log10(Sxw),'linewidth',2); xlim([0 15000]); ylim([-6 1.5]);
xlabel('Frequency (kHz)'); ylabel('Estimated spectrum'); yticks([-6 -4 -2 0]); yticklabels({'10^{-6}' '10^{-4}' '10^{-2}' '10^{0}'}); xticklabels([0 5 10 15]); grid on;
set(gca,'fontsize',18);


% Wavelet transforms:
fmax = 5000;
fmin = 150;
Ms = 250;
freqdisp = [5 4 3 2 1 0.5 0.2]*1e3; % Displayed frequencies in kHz

figure; colormap(flipud(gray));
subplot('Position', [0.05 0.07, 0.45, 0.7]);
display_cwt_JEFAS(y,Fs,fmin,fmax,Ms,wav_typ,wav_param,freqdisp);
xlabel('Time (s)'); ylabel('Frequency (kHz)');
yticklabels(freqdisp/1e3);
title('Original signal');
set(gca,'fontsize',18);

subplot('Position', [0.54 0.07, 0.45, 0.7]);
display_cwt_JEFAS(y,Fs,fmin,fmax,Ms,wav_typ,wav_param,freqdisp);
xlabel('Time (s)'); %ylabel('Frequency (kHz)');
yticklabels(freqdisp/1e3);
title('Stationarized signal');
set(gca,'fontsize',18);
