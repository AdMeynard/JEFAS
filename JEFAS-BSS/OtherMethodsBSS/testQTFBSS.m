clear all; close all;
load('../../signals/synthetic_NonstatMixture.mat');
% addpath(genpath('tftb'));

[N,T] = size(z);

pp = 10; % subsampling for QTF
tt = 1:pp:T;
TT = length(tt);

%% Modified Wigner distribution
tfr1 = tfrpwv([hilbert(y(1,:)).', hilbert(y(1,:)).'],tt,TT);
tfr2 = tfrpwv([hilbert(y(2,:)).', hilbert(y(2,:)).'],tt,TT);
[tfr12,t,f] = tfrpwv([hilbert(y(1,:)).', hilbert(y(2,:)).'],tt,TT);
figure; title('Modified Wigner distribution of the sources');
subplot(2,2,1); imagesc(t,f,abs(tfr1)); axis xy; title('sig 1');
subplot(2,2,2); imagesc(t,f,abs(tfr12)); axis xy; title('sig 1 vs sig2');
subplot(2,2,3); imagesc(t,f,abs(tfr12)); axis xy; title('sig 1 vs sig2');
subplot(2,2,4); imagesc(t,f,abs(tfr2)); axis xy; title('sig 2');

%% Selection of the points
eps1 = 10;
eps2 = 0.1;
eps3 = 0.01;
eps4 = 500;
[ZDmat, Dmat, ptszd, ptsd] = selecpts(z,pp,eps1,eps2, eps3, eps4); % select_pts

Mzd = zeros(TT);
Mzd(ptszd) = 1;
Md = zeros(TT);
Md(ptsd) = 1;
figure; 
subplot(1,2,1); imagesc(t,f,abs(Mzd)); axis xy; title('Zero-diagonalization');
subplot(1,2,2); imagesc(t,f,abs(Md)); axis xy; title('Diagonalization');

%% Algebraic criterion to estimate the mixing matrix 

nn = 2;
Aest = BSS_QTF(z, pp, eps3, eps4, nn); 
display(Aest);

yest = Aest\z;
tfrz1 = tfrpwv([hilbert(z(1,:)).', hilbert(z(1,:)).'],tt,TT);
tfrz2 = tfrpwv([hilbert(z(2,:)).', hilbert(z(2,:)).'],tt,TT);
tfryh1 = tfrpwv([hilbert(yest(1,:)).', hilbert(yest(1,:)).'],tt,TT);
tfryh2 = tfrpwv([hilbert(yest(2,:)).', hilbert(yest(2,:)).'],tt,TT);
figure; title('Modified Wigner distribution of the observationq and estimated sources');
subplot(2,2,1); imagesc(t,f,abs(tfrz1)); axis xy; title('obs. #1');
subplot(2,2,2); imagesc(t,f,abs(tfrz2)); axis xy; title('obs. #2');
subplot(2,2,3); imagesc(t,f,abs(tfryh1)); axis xy; title('est. source  #1');
subplot(2,2,4); imagesc(t,f,abs(tfryh2)); axis xy; title('est. source #2');

