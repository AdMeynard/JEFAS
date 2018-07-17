%% Spectral estimation from a wind sound recording
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
addpath('cwt');
addpath('deform_estimation');
addpath('analysis');

load('signals/wind');
T = length(y);
%% Joint estimation

Dt = 400; % temporal subsampling for the deformation estimation
dgamma0 = ones(1,T); % gamma'(t) initialization
a0 = ones(1,T); % a(t) initialization

wav_typ = 'sharp'; % wavelet type (cf. cwt.m)
wav_paramWP = 20; % corresponding parameter for warping estimation
wav_param = 500; % corresponding parameter for spectrum and AM estimations

NbScales = 125;
scalesAM = 2.^(linspace(2.5,6,NbScales));
subrate = 3; % subsampling step for the scales to ensure the covariance invertibility
scalesWP = scalesAM(1:subrate:end);

r = 1e-5; % regularization parameter

stopWP = 2e-2; % minimal gap between two steps in the gradient
itWP = 6; % number of gradient iterations

Nf = 2500; % number of frequencies for spectrum estimation
NbScalesS = 110;
scalesS = 2.^(linspace(-1,7,NbScalesS)); % for spectrum estimation

Nit = 10; % maximal number of iterations in the joint estimation
stop_crit = 25e-3; % relative update threshold

paramWAV = {wav_typ,wav_param,wav_paramWP};
paramWP = {scalesWP,itWP,stopWP};
paramS = {scalesS,Nf};

paramAM = {'AM',scalesAM,r}; % model with time warping and amplitude modulation
tic;
[aML, dgammaML, Sx, evol_crit] = estim_altern(y,Dt,dgamma0,a0,paramWAV,paramWP,paramAM,paramS,stop_crit,Nit);
toc;


%% Analysis
t = 0:(1/Fs):((T-1)/Fs);
figure;
subplot(2,1,1);plot(t,dgammaML,'linewidth',2); 
ylabel('Estimated log(\gamma''(t))'); axis tight; grid on; ylim([0.5 1.5]);
%set(gca,'FontSize',24);
subplot(2,1,2);plot(t,aML,'linewidth',2); 
xlabel('Time (s)'); ylabel('Estimated a^2(t)'); axis tight; grid on; ylim([0 2]);
%set(gca,'FontSize',24);

z = statAMWP(y,aML,dgammaML);

alpha = 15;
Nff = 50000;
Sxw = estim_spec(z,Nff,alpha);
freq = linspace(0,Fs,Nff);

figure;
semilogy(freq,Sxw,'linewidth',2); 
xlabel('Frequency (Hz)'); ylabel('Estimated spectrum'); axis tight; grid on;
xlim([0 3000]); 
%set(gca,'FontSize',24);

scalesdisp = 2.^(linspace(0.5,3.3,250));
dt = 5;
Wy = cwt(y(1:dt:end),scalesdisp,wav_typ,wav_param);
Wz = cwt(z(1:dt:end),scalesdisp,wav_typ,wav_param);
figure;
subplot(1,2,1); imagesc(abs(Wy));
subplot(1,2,2); imagesc(abs(Wz));

figure;
imagesc(t(1:dt:end),log2(scalesdisp),abs(Wz));
nu0 = Fs/4/dt;
sobs = cellfun(@str2num,get(gca,'yticklabel'));
fobs = round(nu0./2.^sobs);
set(gca,'yticklabel',fobs); %set(gca,'FontSize',26);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); colormap(flipud(gray));