function [c, ceq] = norm_row(B)
%NORM_ROW  Constraint of normalization on the row of a matrix for fmincon
% usage:	[c, ceq] = norm_row(B)
% 
% Input:
%   B: matrix to be normalized
% 
% Output:
%   c: void (required for fmincon)
%   ceq: difference between the norm of the rows of B and one
%
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
% Created: 2018-10-30

P = size(B,1);

v = sqrt(sum(abs(B).^2,2));
N = length(v);
ceq = v - ones(P,1);
c=[];