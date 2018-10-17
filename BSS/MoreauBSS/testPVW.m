clear all; close all;
load('../sig_compBSS.mat');
addpath(genpath('tftb'));

[N,T] = size(z);

pp = 10; % subsampling for QTF
tt = 1:pp:T;
TT = length(tt);

%% Pseudo Wigner Ville
tfr1 = tfrpwv([hilbert(y(1,:)).', hilbert(y(1,:)).'],tt,TT);
tfr2 = tfrpwv([hilbert(y(2,:)).', hilbert(y(2,:)).'],tt,TT);
[tfr12,t,f] = tfrpwv([hilbert(y(1,:)).', hilbert(y(2,:)).'],tt,TT);
figure; title('Pseudo Wigner Ville des sources');
subplot(2,2,1); imagesc(t,f,abs(tfr1)); axis xy; title('sig 1');
subplot(2,2,2); imagesc(t,f,abs(tfr12)); axis xy; title('sig 1 vs sig2');
subplot(2,2,3); imagesc(t,f,abs(tfr12)); axis xy; title('sig 1 vs sig2');
subplot(2,2,4); imagesc(t,f,abs(tfr2)); axis xy; title('sig 2');

%% Points d'interet
eps1 = 10;
eps2 = 0.1;
eps3 = 0.1;
eps4 = 20;
[ZDmat, Dmat, ptszd, ptsd] = selecpts(z,pp,eps1,eps2, eps3, eps4); % select_pts

nn = 9;
Aest = BSS_Moreau(Dmat,N,nn); % BSS dessus
display(Aest);

Mzd = zeros(TT);
Mzd(ptszd) = 1;
Md = zeros(TT);
Md(ptsd) = 1;
figure; 
subplot(1,2,1); imagesc(t,f,abs(Mzd)); axis xy; title('Zéro diagonalisation');
subplot(1,2,2); imagesc(t,f,abs(Md)); axis xy; title('Diagonalisation');

%% Critère algébrique A_DVS
% Kd = size(Dmat,3);
% for k= 1:Kd
%     Dx = real(Dmat(:,:,k));
%     [U,~,~] = svd(Dx);
%     colA(:,k) = U(:,1);
% end
% figure;plot(colA(1,:),colA(2,:),'+')
% [idx,C] = kmeans(colA',nn);
% scatter(colA(1,:),colA(2,:),10,idx)
%  
% nC = histc(idx(:),1:nn);
% for k=1:2
%     A(:,k) = C(nC==max(nC),:).';
%     nC(nC==max(nC)) = 0;
% end
% display(A);

yest = Aest\z;
tfrz1 = tfrpwv([hilbert(z(1,:)).', hilbert(z(1,:)).'],tt,TT);
tfrz2 = tfrpwv([hilbert(z(2,:)).', hilbert(z(2,:)).'],tt,TT);
tfryh1 = tfrpwv([hilbert(yest(1,:)).', hilbert(yest(1,:)).'],tt,TT);
tfryh2 = tfrpwv([hilbert(yest(2,:)).', hilbert(yest(2,:)).'],tt,TT);
figure; title('Pseudo Wigner Ville des observation et sources estimées');
subplot(2,2,1); imagesc(t,f,abs(tfrz1)); axis xy; title('obs 1');
subplot(2,2,2); imagesc(t,f,abs(tfrz2)); axis xy; title('obs 2');
subplot(2,2,3); imagesc(t,f,abs(tfryh1)); axis xy; title('source est. 1');
subplot(2,2,4); imagesc(t,f,abs(tfryh2)); axis xy; title('source est. 2');

