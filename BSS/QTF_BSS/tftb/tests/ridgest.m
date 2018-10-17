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
%function ridgest
%RIDGEST Unit test for the function RIDGES.

%	O. Lemoine - April 1996.

N=128; 
t=1:N;

% time localization
t0=60; sig=((1:N)'==t0);
[tfr,rtfr,hat]=tfrrsp(sig);
[ptt,ptf]=ridges(tfr,hat,t,'tfrrsp'); 
if any(ptt~=t0)|any(abs(ptf'-(1:N)/N)>sqrt(eps)),
 error('ridges test 1 failed');
end;


% frequency localization
f0=30;
sig=fmconst(N,f0/128);
[tfr,rtfr,hat]=tfrrsp(sig);
[ptt,ptf]=ridges(tfr,hat,t,'tfrrsp'); 
if any((ptf(3:106)-(f0+1)/N)>sqrt(eps)),
 error('ridges test 2 failed');
end;
if any(abs(ptt(3:106)-(13:116)')>sqrt(eps)),
 error('ridges test 3 failed');
end;


N=117; 
t=1:N;

% time localization
t0=53; sig=((1:N)'==t0);
[tfr,rtfr,hat]=tfrrsp(sig);
[ptt,ptf]=ridges(tfr,hat,t,'tfrrsp'); 
if any(ptt~=t0)|any(abs(ptf'-(1:N)/N)>sqrt(eps)),
 error('ridges test 4 failed');
end;


% frequency localization
f0=31;
sig=fmconst(N,f0/N);
[tfr,rtfr,hat]=tfrrsp(sig);
[ptt,ptf]=ridges(tfr,hat,t,'tfrrsp'); 
if any((ptf(3:95)-(f0+1)/N)>sqrt(eps)),
 error('ridges test 5 failed');
end;
if any(abs(ptt(3:95)-(11:103)')>sqrt(eps)),
 error('ridges test 6 failed');
end;

