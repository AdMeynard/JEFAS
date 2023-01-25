function thetaAM0 = baselineAMest(Wy)
%BASELINEAMEST	Compute the local scale estimator of a^2
% usage:	thetaAM0 = baselineAMest(Wy)
%
% Input:
%   Wy: wavelet transform of the deformed signal
% 
% Output:
%   thetaAM0 : local scale estimator of a^2

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


scalogram = abs(Wy).^2;
est = sum(scalogram);
thetaAM0 = est./mean(est);  % such as <a^2>=1

end