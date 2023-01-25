function [ordered_newsig, stablematch] = reordersig(oldsig,newsig)
% REORDERSIG Prevent estimated sources order from being changed througth JEFAS-BSS iterations
% usage:	[ordered_newsig, stablematch] = reordersig(oldsig,newsig)
%
% Input:
%   oldsig: previous estimated sources
%   newsig: current estimated sources
%
% Output:
%   ordered_newsig: reordered current estimated sources
%   stablematch: index ordering the estimated sources

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

N = size(newsig,1);

for n=1:N
    x = newsig(n,:);
    for k=1:N
        y = oldsig(k,:);
        M(n,k) = max(abs(xcorr(x,y)));
    end
end

[~,Prefm] = sort(-M,2);
[~,Prefw] = sort(-M,1);
Prefw = Prefw.';

stablematch = galeshapley(N,Prefm,Prefw);
ordered_newsig = newsig(stablematch,:);
