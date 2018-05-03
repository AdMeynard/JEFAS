function W = cwt(sig,scales,wav_typ,wav_par)
% CWT:	Continuous wavelet transform with analytic derivative of Gaussian
%       wavelet or the sharp wavelet
% usage:	W = cwt(sig,scales,wav_typ,k)
%
% Input:
%   sig: vector containing input signal
%   scales : vector of scales
%   wav_typ: wavelet type:
%     wav_typ=0: sharp
%     wav_typ=1: derivative of gaussian
%   wav_par: parameter depending on the wavelet type
%     wav_typ=0: wav_par = -ln(epsilon)>0 where epsilon=value of \hat{\psi}(Fs/2)
%     wav_typ=1: wav_par = number of vanishing moments
%
% Output:
%   W: continuous wavelet transform (complex valued)

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
% Created: 2017-12-19

siglength = length(sig);
fff = (0:(siglength-1))*2*pi/siglength;

fsig = fft(sig(:)).';

scales = scales(:);

tmp = scales * fff;

% Generate wavelets in the Fourier domain
switch wav_typ
    case 0
        fpsi = exp(-2*wav_par*((pi./(2*tmp)) + (2*tmp/pi) - 2));

    case 1
        Cst = 4*wav_par/(pi^2); % such as argmax fpsi = pi/2
        K = (2/pi)^wav_par*exp(wav_par/2); % such as fpsi(pi/2) = 1
        t_powV = tmp.^wav_par;
        fpsi = K*t_powV.* exp(-Cst*tmp.^2/2);
    
    otherwise
        error('Unexpected wavelet type. CWT not computed')    
end

U = bsxfun(@times,sqrt(scales),fpsi);
fTrans = bsxfun(@times,U,fsig);
W = ifft(fTrans,[],2);
    
end