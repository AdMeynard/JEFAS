function [W,t,s] = display_cwt_JEFAS(sig,Fs,fmin,fmax,Ms,wav_typ,wav_par,varargin)
% display_cwt_JEFAS:	Display the continuous wavelet transform with analytic derivative of Gaussian
%       wavelet or the sharp wavelet
% usage:	W = display_cwt_JEFAS(sig,Fs,fmin,fmax,Ms,wav_typ,wav_par,freqtick)
%   /!\ If you do not want to display the CWT, replace freqtick by a 0 argument
%
% Input:
%   sig: vector containing the input signal
%   sig: sampling frequency
%   fmin: minimum frequency to display
%   fmax: maximum frequency to display
%   Ms: length of the scales vector (= number of rows of the CWT) 
%   wav_typ: wavelet type:
%     wav_typ='sharp': The sharp wavelet
%     wav_typ='dgauss': analytic derivative of gaussian
%   wav_par: parameter depending on the wavelet type
%     if wav_typ='sharp': wav_par = -ln(epsilon)>0 where epsilon=value of \hat{\psi}(Fs/2)
%     if wav_typ='dgauss': wav_par = number of vanishing moments
%   freqtick (optional): frequency to tick on the y-axis
%
% Output:
%   W: wavelet transform coefficients matrix (complex valued)

% Copyright (C) 2017 Adrien MEYNARD
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
% Created: 2019-03-13

N = length(sig);
t = linspace(0,(N-1)/Fs,N);

% Scales:
xi0 = Fs/4; % wavelet central frequency
smin = log2( xi0/fmax );
smax = log2( xi0/fmin );
scales = 2.^(linspace(smin,smax,Ms));
s = log2(scales);

% Wavelet transform:
W = cwt_JEFAS(sig,scales,wav_typ,wav_par);


% Display:
if ( nargin==8 && isequal(varargin{1},0) )
    return
else
    if nargin==8
        freqtick = sort(varargin{1},'descend');
    else
        freqtick = floor(fmax*(fmin/fmax).^(linspace(0,1,5)));
    end
    sdisp = log2(xi0./freqtick); % corresponding log-scales

    imagesc(t,s,log1p(abs(W))/0.1);
    yticks(sdisp); yticklabels(freqtick);
end

end