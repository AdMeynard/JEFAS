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
%function noisecut
%NOISECUT Unit test for the function NOISECU.

%	O. Lemoine - May 1996.

N=32768;
sig=noisecu(N);

% Mean
Mean=mean(sig);
if abs(Mean)>10/sqrt(N),
 error('noisecu test 1 failed');
end


% Variance
Var=std(sig).^2;
if abs(Var-1)>10/sqrt(N),
 error('noisecu test 2 failed');
end


% histogram
sig=noisecu(N);
Nh=100;
h=hist(real(sig),Nh); h=h/mean(h);
pdf=ones(1,Nh);
if any(abs(h-pdf).^2>10/sqrt(N)),
 error('noisecu test 3 failed');
end


% For N=1 
N=1; Np=10000;
for k=1:Np,
 sig(k)=noisecu(N);
end
Mean=mean(sig);
if abs(Mean)>10/sqrt(Np),
 error('noisecu test 4 failed');
end
Var=std(sig).^2;
if abs(Var-1)>10/sqrt(Np),
 error('noisecu test 5 failed');
end
h=hist(real(sig),Nh); h=h/mean(h);
if any(abs(h-pdf).^2>10/sqrt(Np)),
 error('noisecu test 6 failed');
end


% For N=2
N=2;
for k=1:2:(Np-1),
 noise=noisecu(N);
 sig(k)=noise(1);
 sig(k+1)=noise(2);
end
Mean=mean(sig);
if abs(Mean)>10/sqrt(Np),
 error('noisecu test 7 failed');
end
Var=std(sig).^2;
if abs(Var-1)>10/sqrt(Np),
 error('noisecu test 8 failed');
end
h=hist(real(sig),Nh); h=h/mean(h);
if any(abs(h-pdf).^2>10/sqrt(Np)),
 error('noisecu test 9 failed');
end

