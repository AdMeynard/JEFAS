function act = sigDetection(y,ratio,scales,wav_typ,wav_param)
% SIGDETECTION Detect the activity of the signal
% usage:	act = sigDetection(y,ratio,scales,wav_typ,wav_param)
%
% Input:
%   y: signal to analyze
%   ratio: threshold do detect activity
%   scales: vector of scales for wavelet transform
%   wav_typ: cf. cwt_JEFAS
%   wav_param: cf. cwt_JEFAS
%
% Output:
%   act: is one when the signal is active, zero elsewhere

% Copyright (C) 2019 Adrien MEYNARD
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
% Created: 2019-03-11

T = length(y);
act = ones(1,T);

Wy = cwt_JEFAS(y,scales,wav_typ,wav_param);
evolPower = sum( abs(Wy).^2 );

medPower = median(evolPower);
threshold = ratio*medPower;

act(evolPower <= threshold) = 0;