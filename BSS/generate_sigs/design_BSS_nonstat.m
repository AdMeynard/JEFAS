addpath('../../Sig_design');
addpath('../../cwt');

clear all; close all;
%% Génération de signaux stationaires
N = 2;
T = 44100;
x = randn(T,N);

[x1,Sx1] = BandPassApprox(x(:,1),floor(0.06*T),floor(0.1*T));
[x2,Sx2] = BandPassApprox(x(:,1),floor(0.15*T),floor(0.2*T));

Sx = [Sx1'; Sx2'];

%% Déformation de ceux-ci

[y1,~,dgamma1] = chirpwarp(x1,1);
[y2,~,dgamma2] = sinewarp(x2,2,0.02);

y = [y1';y2'];
dgamma = [dgamma1'; dgamma2'];

%% Mixing matrix

typmel = 'stat';

switch typmel
    case 'stat'
        A = [1 0.75; -0.5 1];
        z = A*y;
%         save('sig_compBSS','y','z','A','Sx','dgamma');
      
    case 'nonstat'
        pileA = zeros(N,N,T);
        pileB = zeros(N,N,T);
        z = zeros(N,T);
        for t=1:T
            A = [1+0.3*cos(5*pi*t/T)/4 0.75+0.4*cos(3*pi*t/T)/4;
                -0.5+0.5*cos(11*pi*t/T)/4 1+2*0.1*cos(8*pi*t/T)/4];
            B = inv(A);
            z(:,t) = A*y(:,t);
            c(t)=cond(A);
            pileA(:,:,t) = A;
            pileB(:,:,t) = B;
        end
        figure;
        plot(c);
        title('condtionnement de A');
        save('sig_compBSS_nonstat2','y','z','pileA','pileB','Sx','dgamma');
        figure;
        b11=pileB(1,1,:); b11 = b11(:);
        b12=pileB(1,2,:); b12 = b12(:);
        b21=pileB(2,1,:); b21 = b21(:);
        b22=pileB(2,2,:); b22 = b22(:);
        plot(1:T,b11,'b',1:T,b12,'k',1:T,b21,'r',1:T,b22,'g');
        title('éléments de B');
%         save('sig_compBSS_nonstat','y','z','A','Sx','dgamma');
end

scales = 2.^(linspace(-1,4,100));
Wz1 = cwt(z(1,:),scales,'sharp',800);
figure; imagesc(abs(Wz1));
title('Observations 1');

Wz2 = cwt(z(2,:),scales,'sharp',800);
figure; imagesc(abs(Wz2));
title('Observations 2');

