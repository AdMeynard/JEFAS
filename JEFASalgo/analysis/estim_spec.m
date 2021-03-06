function S = estim_spec(x,Nf,alpha)
% ESTIM_SPEC estimation the spectrum of a stationary signal (Welch method)
% usage:	Sz = estim_spec(x,Nf,alpha)
%
% Input:
%   x: stationary signal
%   Nf: length of the frequency vector on which S is computed (between 0
%   and Fs)
%   alpha(>=1): relative to the analysis window length. The higher alpha,
%   the smoother the estimated spectrum  
% 
% Output:
%   S : estimated spectrum of x

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
% Created: 2017-12-20

if alpha<1
    error('alpha must be greater than one')
end

T = length(x);
nwin = floor(T/alpha);
S = pwelch(x,nwin,[],Nf);
S = [S; flipud(S(1:(end-1)))]; % real signal
S = pi*S(1:(end-1))';

end