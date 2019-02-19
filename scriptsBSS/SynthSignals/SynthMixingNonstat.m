addpath('../../cwt');
clear all; close all;

K = 20;
for k= 1:K
%% Stationary sources 
    N = 2;
    T = 44100;
    x = randn(T,N);

    [x1,Sx1] = BandPassApprox(x(:,1),floor(0.03*T),floor(0.05*T));
    [x2,Sx2] = BandPassApprox(x(:,2),floor(0.075*T),floor(0.10*T));

    Sx = [Sx1'; Sx2'];

    %% Time warpings

    [y1,~,dgamma1] = chirpwarp(x1,1);
    [y2,~,dgamma2] = sinewarp(x2,2,0.02);

    sigmay = 5e-2;
    y = [y1';y2'] + sigmay*randn(N,T);
    dgamma = [dgamma1'; dgamma2'];

    %% Mixing matrix

    typmel = 'nonstat';

    switch typmel
        case 'stat'
            A = [1 0.75; -0.5 1];
            z = A*y;

        case 'nonstat'
            heapA = zeros(N,N,T);
            heapB = zeros(N,N,T);
            z = zeros(N,T);
            for t=1:T
                A = [1+0.3*cos(5*pi*t/T) 0.75+0.4*cos(3*pi*t/T);
                    -0.5+0.5*cos(11*pi*t/T) 1+0.1*cos(8*pi*t/T)];
                B = inv(A);
                z(:,t) = A*y(:,t);
                c(t)=cond(A);
                heapA(:,:,t) = A;
                heapB(:,:,t) = B;
            end
%             figure;
%             plot(c);
%             title('Condition number');
%             figure;
%             b11=heapB(1,1,:); b11 = b11(:);
%             b12=heapB(1,2,:); b12 = b12(:);
%             b21=heapB(2,1,:); b21 = b21(:);
%             b22=heapB(2,2,:); b22 = b22(:);
%             plot(1:T,b11,'b',1:T,b12,'k',1:T,b21,'r',1:T,b22,'g');
%             title('B entries');
            
            heapAK(:,:,:,k) = heapA;
            heapBK (:,:,:,k) = heapB;
    end
    
    yK(:,:,k) = y;
    zK(:,:,k) = z;
    SxK (:,:,k) = Sx;
    dgammaK (:,:,k) = dgamma; 
end

%% Analysis
wav_typ = 'sharp';
wav_param = 500;

NbScales = 100;
scales = 2.^(linspace(0,5,NbScales));

Fs = 44100;
t = linspace(0,(T-1)/Fs,T);
nu0 = Fs/4;
freqdisp = [16 8 4 2 1]; % Displayed frequencies in kHz
sdisp = log2(nu0./(1e3*freqdisp)); % corresponding log-scales

figure; title('Sources and observations');
for n=1:N
    Wy = cwt_JEFAS(y(n,:),scales,wav_typ,wav_param);
    subplot(N,N,n); imagesc(t,log2(scales),abs(Wy)); yticks(sdisp); yticklabels(freqdisp);
    ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
    Wz = cwt_JEFAS(z(n,:),scales,wav_typ,wav_param);
    subplot(N,N,N+n); imagesc(t,log2(scales),abs(Wz)); yticks(sdisp); yticklabels(freqdisp);
    xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
end

