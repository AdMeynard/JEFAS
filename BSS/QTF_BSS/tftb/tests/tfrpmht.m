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
%function tfrpmht
%TFRPMHT Unit test for the time frequency representation TFRPMH.

%       O. Lemoine - March 1996. 

% We test each property of the corresponding TFR :

N=128;

% Covariance by translation in time 
t1=55; t2=70; f=0.3;
sig1=amgauss(N,t1).*fmconst(N,f,t1); 
sig2=amgauss(N,t2).*fmconst(N,f,t2); 
tfr1=tfrpmh(sig1);  
tfr2=tfrpmh(sig2);        
[tr,tc]=size(tfr1);
nu=round(f*(tc-1)*2)+1;
tfr=tfr1-tfr2(:,modulo((1:tc)-t1+t2,tc));
if any(any(abs(tfr)>sqrt(eps))),
 error('tfrpmh test 1 failed');
end


% Reality of the TFR
sig=noisecg(N);
tfr=tfrpmh(sig);
if sum(any(abs(imag(tfr))>sqrt(eps)))~=0,
 error('tfrpmh test 2 failed');
end


% Energy conservation
sig=noisecg(N);
tfr=tfrpmh(sig);
Es=norm(sig)^2;
Etfr=sum(mean(tfr));
if abs(Es-Etfr)>sqrt(eps),
 error('tfrpmh test 3 failed');
end


% Time-marginal
sig=noisecg(N);
tfr=tfrpmh(sig);
ip1=abs(sig).^2;
ip2=mean(tfr)';
if any(abs(ip1-ip2)>sqrt(eps)),
 error('tfrpmh test 4 failed');
end


% Conservation of the time support (wide-sense)
sig=[zeros(N/4,1);noisecg(N/2);zeros(N/4,1)];
tfr=tfrpmh(sig);
if sum(any(abs(tfr(:,1:N/4-1))>sqrt(eps))) | ...
   sum(any(abs(tfr(:,(3*N/4+1):N))>sqrt(eps))),
 error('tfrpmh test 5 failed');
end


% time localization
t0=30; sig=((1:N)'==t0);
tfr=tfrpmh(sig);
[ik,jk]=find(tfr~=0.0);
if any(jk~=t0)|any(ik'-(1:N)),
 error('tfrpmh test 6 failed');
end;


% A PMHD with a constant window is a MHD
sig=fmlin(N).*amgauss(N);
tfr1=tfrpmh(sig,1:N,N,ones(2*N+1,1));
tfr2=tfrmh(sig);
if max(max(abs(tfr1-tfr2)))>1e-5,
 error('tfrpmh test 7 failed');
end


N=131;

% Covariance by translation in time 
t1=55; t2=70; f=0.3;
sig1=amgauss(N,t1).*fmconst(N,f,t1); 
sig2=amgauss(N,t2).*fmconst(N,f,t2); 
tfr1=tfrpmh(sig1);  
tfr2=tfrpmh(sig2);        
[tr,tc]=size(tfr1);
nu=round(f*(tc-1)*2)+1;
tfr=tfr1-tfr2(:,modulo((1:tc)-t1+t2,tc));
if any(any(abs(tfr)>sqrt(eps))),
 error('tfrpmh test 8 failed');
end


% Reality of the TFR
sig=noisecg(N);
tfr=tfrpmh(sig);
if sum(any(abs(imag(tfr))>sqrt(eps)))~=0,
 error('tfrpmh test 9 failed');
end


% Energy conservation
sig=noisecg(N);
tfr=tfrpmh(sig);
Es=norm(sig)^2;
Etfr=sum(mean(tfr));
if abs(Es-Etfr)>sqrt(eps),
 error('tfrpmh test 10 failed');
end


% Time-marginal
sig=noisecg(N);
tfr=tfrpmh(sig);
ip1=abs(sig).^2;
ip2=mean(tfr)';
if any(abs(ip1-ip2)>sqrt(eps)),
 error('tfrpmh test 11 failed');
end


% Conservation of the time support (wide-sense)
sig=[zeros(round(N/4),1);noisecg(round(N/2));zeros(round(N/4),1)];
tfr=tfrpmh(sig);
if sum(any(abs(tfr(:,1:round(N/4)-1))>sqrt(eps))) | ...
   sum(any(abs(tfr(:,(round(3*N/4)+2):N))>sqrt(eps))),
 error('tfrpmh test 12 failed');
end


% time localization
t0=30; sig=((1:N)'==t0);
tfr=tfrpmh(sig);
[ik,jk]=find(tfr~=0.0);
if any(jk~=t0)|any(ik'-(1:N)),
 error('tfrpmh test 13 failed');
end;


% A PMHD with a constant window is a MHD
sig=fmlin(N).*amgauss(N);
tfr1=tfrpmh(sig,1:N,N,ones(2*N+1,1));
tfr2=tfrmh(sig);
if max(max(abs(tfr1-tfr2)))>1e-5,
 error('tfrpmh test 14 failed');
end

