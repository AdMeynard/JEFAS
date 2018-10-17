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
%function integ2dt
%INTEG2DT Unit test for the function INTEG2D.

%	O. Lemoine - March 1996.

N=128;

% Integration of a constant matrix
MAT=ones(N);
E = integ2d(MAT);
if abs(E-(N-1)^2)>sqrt(eps),
 error('integ2d test 1 failed');
end

% Integration of a sinusoidal matrix
x=(0:N);
MAT=sin(2*pi*x/N)'*sin(2*pi*x/N);
E = integ2d(MAT);
if abs(E)>sqrt(eps),
 error('integ2d test 2 failed');
end

% Integration of a scalogram
sig=altes(N,0.1,0.45,10000);
Es=norm(sig); 
[TFR,T,F] = tfrscalo(sig,1:N,0,.01,.49,2*N) ;
E = integ2d(TFR,T,F);
if abs(E-1)>sqrt(eps),
 error('integ2d test 3 failed');
end


N=117;

% Integration of a constant matrix
MAT=ones(N);
E = integ2d(MAT);
if abs(E-(N-1)^2)>sqrt(eps),
 error('integ2d test 4 failed');
end

% Integration of a sinusoidal matrix
x=(0:N);
MAT=sin(2*pi*x/N)'*sin(2*pi*x/N);
E = integ2d(MAT);
if abs(E)>sqrt(eps),
 error('integ2d test 5 failed');
end

% Integration of a scalogram
sig=altes(N,0.1,0.45,10000);
Es=norm(sig); 
[TFR,T,F] = tfrscalo(sig,1:N,0,.01,.49,2*N) ;
E = integ2d(TFR,T,F);
if abs(E-1)>sqrt(eps),
 error('integ2d test 6 failed');
end

