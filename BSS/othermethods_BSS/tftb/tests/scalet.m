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
%function scalet
%SCALET	Unit test for the function SCALE.

%	O. Lemoine - March 1996.

N=128;

t=N/2; f=0.4; T=2*sqrt(N); a=2;
sig1=amgauss(N,t,T).*fmconst(N,f,t);
sig2=amgauss(a*N,a*t,T*a).*fmconst(a*N,f/a,a*t)/sqrt(a);
sig1sc=scale(sig1,2,0.01,.49,2*N);
Diff=abs(sig1sc-sig2);
if any(Diff>1e-7),
  error('scale test 1 failed');
end

 
t=round(N/3); f=0.3; T=sqrt(N); a=3;
sig1=amgauss(N,t,T).*fmconst(N,f,t);
sig2=amgauss(a*N,N,T*a).*fmconst(a*N,f/a,N)/sqrt(a);
sig1sc=scale(sig1,a,0.01,.49,2*N);
Diff=sig1sc-sig2;
if any(abs(Diff)>1e-7),
  error('scale test 2 failed');
end


t=N/2; f=0.15; T=2*sqrt(N); a=1/2;
sig1=amgauss(N,t,T).*fmconst(N,f,t);
sig2=sig1(1:1/a:N)/sqrt(a);
sig1sc=scale(sig1,a,0.01,.49,2*N);
Diff=sig1sc(2:N/2)-sig2(1:N/2-1);
if any(abs(Diff)>1e-7),
  error('scale test 3 failed');
end


N=125;

t=round(N/2); f=0.4; T=2*sqrt(N); a=2;
sig1=amgauss(N,t,T).*fmconst(N,f,t);
sig2=amgauss(a*N,a*t,T*a).*fmconst(a*N,f/a,a*t)/sqrt(a);
sig1sc=scale(sig1,2,0.01,.49,2*N);
Diff=abs(sig1sc(1:N-1)-sig2(2:N));
if any(Diff>1e-7),
  error('scale test 4 failed');
end

 
t=round(N/3); f=0.3; T=sqrt(N); a=3;
sig1=amgauss(N,t,T).*fmconst(N,f,t);
sig2=amgauss(a*N,N,T*a).*fmconst(a*N,f/a,N)/sqrt(a);
sig1sc=scale(sig1,a,0.01,.49,2*N);
Diff=sig1sc(1:a*N-1)-sig2(2:a*N);
if any(abs(Diff)>1e-7),
  error('scale test 5 failed');
end

