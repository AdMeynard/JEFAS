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
%function tfrideat
%TFRIDEAT Unit test for the function TFRIDEAL.

%	O. Lemoine - March 1996.

N=128;

[sig1,ifl1]=fmlin(N); [sig2,ifl2]=fmsin(N);
[sig3,ifl3]=fmhyp(N,[1 0.5],[N/2 0.1]); 
tfr=tfrideal([ifl1;ifl2;ifl3]); 
if any(any(tfr~=1 & tfr~=0)),
 error('tfrideal test 1 failed');
elseif any(sum(tfr)~=1),
 error('tfrideal test 2 failed');
end

iflaws1=[ifl1;ifl2;ifl3]';
for t=1:3*N,
  [Max iflaws2(t)]=max(tfr(:,t));
end
iflaws2=(iflaws2-1)*0.5/(3*N-1);

if any(abs(iflaws1-iflaws2)>1/(3*N)), 
 error('tfrideal test 3 failed');
end


clear; N=111;

[sig1,ifl1]=fmlin(N); [sig2,ifl2]=fmsin(N);
[sig3,ifl3]=fmhyp(N,[1 0.5],[round(N/2) 0.1]); 
tfr=tfrideal([ifl1;ifl2;ifl3]); 
if any(any(tfr~=1 & tfr~=0)),
 error('tfrideal test 4 failed');
elseif any(sum(tfr)~=1),
 error('tfrideal test 5 failed');
end

iflaws1=[ifl1;ifl2;ifl3]';
for t=1:3*N,
  [Max iflaws2(t)]=max(tfr(:,t));
end
iflaws2=(iflaws2-1)*0.5/(3*N-1);

if any(abs(iflaws1-iflaws2)>1/(3*N)), 
 error('tfrideal test 6 failed');
end
