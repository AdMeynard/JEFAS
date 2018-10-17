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
%function tfrrpmht
%TFRRPMHT Unit test for the function TFRRPMH.

%	O. Lemoine - April 1996. 


N=128;

% Reality of the TFR
sig=noisecg(N);
[tfr,rtfr]=tfrrpmh(sig);
if sum(any(abs(imag(rtfr))>sqrt(eps)))~=0,
 error('tfrrpmh test 1 failed');
end


% Energy conservation
sig=fmlin(N);
[tfr,rtfr]=tfrrpmh(sig);
Es=norm(sig)^2;
Etfr=sum(mean(rtfr));
if abs(Es-Etfr)>sqrt(eps),
 error('tfrrpmh test 2 failed');
end


% time localization
t0=30; sig=((1:N)'==t0);
[tfr,rtfr]=tfrrpmh(sig);
[ik,jk]=find(abs(rtfr)>sqrt(eps));
if any(jk~=t0)|any(ik'-(1:N)),
 error('tfrrpmh test 3 failed');
end;


% frequency localization
f0=30;
sig=fmconst(N+6,f0/N);
[tfr rtfr]=tfrrpmh(sig,N/2+2,N,tftb_window(N+1,'rect'));
if any(find(rtfr>2*max(rtfr)/N)~=f0+1)|(abs(mean(rtfr)-1.0)>sqrt(eps)),
 error('tfrrpmh test 4 failed');
end;


N=127;

% Reality of the TFR
sig=noisecg(N);
[tfr,rtfr]=tfrrpmh(sig);
if sum(any(abs(imag(rtfr))>sqrt(eps)))~=0,
 error('tfrrpmh test 5 failed');
end


% Energy conservation
sig=fmlin(N);
[tfr,rtfr]=tfrrpmh(sig);
Es=norm(sig)^2;
Etfr=sum(mean(rtfr));
if abs(Es-Etfr)>sqrt(eps),
 error('tfrrpmh test 6 failed');
end


% time localization
t0=30; sig=((1:N)'==t0);
[tfr,rtfr]=tfrrpmh(sig);
[ik,jk]=find(abs(rtfr)>sqrt(eps));
if any(jk~=t0)|any(ik'-(1:N)),
 error('tfrrpmh test 7 failed');
end;


% frequency localization
f0=30;
sig=fmconst(N+6,f0/N);
[tfr rtfr]=tfrrpmh(sig,round(N/2)+2,N,tftb_window(N,'rect'));
if any(find(rtfr>2*max(rtfr)/N)~=f0+1)|(abs(mean(rtfr)-1.0)>sqrt(eps)),
 error('tfrrpmh test 8 failed');
end;

