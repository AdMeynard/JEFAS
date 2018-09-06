function y = icwt_JEFAS(W,scales,wav_typ,wav_par)
%ICWT:	Continuous inverse wavelet transform with analytic derivative of
%       Gaussian or the sharp wavelet
% usage:	y = icwt_JEFAS(W,scales,wav_typ,wav_par)
%
%
% Input:
%   W: wavelet transform coefficients matrix
%   scales : vector of scales
%   wav_typ: wavelet type:
%     wav_typ='sharp': The sharp wavelet
%     wav_typ='dgauss': analytic derivative of gaussian
%   wav_par: parameter depending on the wavelet type
%     if wav_typ='sharp': wav_par = -ln(epsilon)>0 where epsilon=value of \hat{\psi}(Fs/2)
%     if wav_typ='dgaudd: wav_par = value of \hat{\psi}(Fs/2)
%
% Output:
%   y: reconstructed signal

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

siglength = size(W,2);
fff = (0:(siglength-1))*2*pi/siglength;

scales = scales(:);
ds = log(scales(2)) - log(scales(1));

tmp = scales * fff;

% Generate wavelets in the Fourier domain
switch wav_typ
    case 'sharp'
        fpsi = exp(-2*wav_par*((pi./(2*tmp)) + (2*tmp/pi) - 2));
        Cpsi = 0.88./wav_par.^0.5; % empirical

    case 'dgauss'
        Cst = 4*wav_par/(pi^2); % such that argmax fpsi = pi/2
        K = (2/pi)^wav_par*exp(wav_par/2); % such that fpsi(pi/2) = 1
        t_powV = tmp.^wav_par;
        fpsi = K*t_powV.* exp(-Cst*tmp.^2/2);
        Cpsi = (K^2/(2*Cst^wav_par))*gamma(wav_par); % theoretical
    
    otherwise
        error('Unexpected wavelet type. ICWT not computed')    
end

fy = fpsi .* fft(W,[],2);
fy = 2*real(ifft(fy,[],2));
fy = bsxfun(@times,1./scales.^0.5,fy);
y = (1/Cpsi)*sum(fy)*ds;
    
end