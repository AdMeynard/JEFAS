%% Spectral estimation from a synthesized rapidly wind sound
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
% Created: 2022-03-07

%% Load signal

clear all; close all; clc;

warning off;
addpath('../cwt');
addpath(genpath('../JEFASalgo'));

load('../signals/wind4JEFASS.mat');
T = length(y);
%% Joint estimation

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

tic;
[aML, dgammaML, SxML, evol_crit] = estim_altern(y,Dt,ratio,dgamma0,a0,paramWAV,paramWP,paramAM,paramS,stop_crit,Nit);
toc;

save('results/JEFAS_wind4JEFASS','aML', 'dgammaML','SxML','y','Fs')


%% Analysis

% Deformations:
t = 0:(1/Fs):((T-1)/Fs);
figure;
plot(t,dgamma,'k--',t,dgammaML,'r','linewidth',2);
xlabel('Time (s)'); ylabel("\gamma'");
legend('Ground truth function','JEFAS estimate');
axis tight; grid on; ylim([0.5 1.5]);
set(gca,'fontsize',18);
