%% Cross-synthetis of sounds
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
addpath(genpath('../JEFASalgo'));

load('results/res_wind')
awind = aML;
dgammawind = dgammaML;
zwind = statAMWP(y,aML,dgammaML); % stationarization

load('results/res_sing')
asing = aML;
dgammasing = dgammaML;
zsing = statAMWP(y,aML,dgammaML);

load('results/res_f1')
af1 = aML;
dgammaf1 = dgammaML;
zf1 = statAMWP(y,aML,dgammaML);

%% Synthesis

ywind2sing = deformAMWP(zwind,asing,dgammasing);
ywind2f1 = deformAMWP(zwind,af1,dgammaf1);
ysing2wind = deformAMWP(zsing,awind,dgammawind);
ysing2f1 = deformAMWP(zsing,af1,dgammaf1);
yf12wind = deformAMWP(zf1,awind,dgammawind);
yf12sing = deformAMWP(zf1,asing,dgammasing);