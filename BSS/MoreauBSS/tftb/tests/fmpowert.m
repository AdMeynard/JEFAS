%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%function fmpowert
%FMPOWERT Unit test for the function FMPOWER.

%	O. Lemoine - February, April 1996.

N=256;
t=1:N;

% Mentionning F0 and C
K=1.1; F0=0.05; C=.3;
[x,iflaw]=fmpower(N,K,[F0,C]);
if any(abs(iflaw-(F0+C./t.^K)')>sqrt(eps))~=0,
  error('fmpower test 1 failed');
end

% Mentionning P1 and P2
K=1.8; P1=[1,.1]; P2=[200,.4];
[x,iflaw]=fmpower(N,K,P1,P2);
if abs(iflaw(P1(1))-P1(2))>sqrt(eps) | abs(iflaw(P2(1))-P2(2))>sqrt(eps),
  error('fmpower test 2 failed');
end

% Other values of K, F0 and C
K=-2.1; F0=0.01; C=3.5e-6;
[x,iflaw]=fmpower(N,K,[F0,C]);
if any(abs(iflaw-(F0+C./t.^K)')>sqrt(eps))~=0,
  error('fmpower test 3 failed');
end

% Other values of K, P1 and P2
K=-1.3; P1=[1,.42]; P2=[N,0.03];
[x,iflaw]=fmpower(N,K,P1,P2);
if abs(iflaw(P1(1))-P1(2))>sqrt(eps) | abs(iflaw(P2(1))-P2(2))>sqrt(eps),
  error('fmpower test 4 failed');
end


N=251;
t=1:N;

% Mentionning F0 and C
K=1.1; F0=0.05; C=.3;
[x,iflaw]=fmpower(N,K,[F0,C]);
if any(abs(iflaw-(F0+C./t.^K)')>sqrt(eps))~=0,
  error('fmpower test 5 failed');
end

% Mentionning P1 and P2
K=1.8; P1=[1,.1]; P2=[200,.4];
[x,iflaw]=fmpower(N,K,P1,P2);
if abs(iflaw(P1(1))-P1(2))>sqrt(eps) | abs(iflaw(P2(1))-P2(2))>sqrt(eps),
  error('fmpower test 6 failed');
end

% Other values of K, F0 and C
K=-2.1; F0=0.01; C=3.5e-6;
[x,iflaw]=fmpower(N,K,[F0,C]);
if any(abs(iflaw-(F0+C./t.^K)')>sqrt(eps))~=0,
  error('fmpower test 7 failed');
end

% Other values of K, P1 and P2
K=-1.3; P1=[1,.42]; P2=[N,0.03];
[x,iflaw]=fmpower(N,K,P1,P2);
if abs(iflaw(P1(1))-P1(2))>sqrt(eps) | abs(iflaw(P2(1))-P2(2))>sqrt(eps),
  error('fmpower test 8 failed');
end

