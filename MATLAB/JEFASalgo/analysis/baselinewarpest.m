function thetaWP0 = baselinewarpest(Wy,scales)
%BASELINEWARPEST	Compute the local scale estimator of log2(gamma')
% usage:	thetaWP0 = baselinewarpest(Wy,scales)
%
% Input:
%   Wy: wavelet transform of the deformed signal
%   scales: scales on which the CWT is calculated
% 
% Output:
%   thetaWP0 : local scale estimator of log2(gamma')

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

scales_log = log2(scales);

scalogram = abs(Wy).^2;
smoy = scales_log(:)' * scalogram;
normwx = sum(scalogram);
est = smoy ./normwx;

thetaWP0 = - est - log2(mean(2.^(-est))); % such that <gamma'>=1
end