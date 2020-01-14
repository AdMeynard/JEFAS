addpath('../../cwt');
addpath('../../signals');
clear all; close all;

load wind
[P,Q] = rat(44100/Fs);
y1 = resample(y,P,Q).';
load sing
y2=y(:).';

T = min(length(y1),length(y2));
y1 = y1(1:T)/std(y1(1:T)) ;
y2 = y2(1:T)/std(y2(1:T)) ;

y = [y1; y2];
[N,T] = size(y);

f0 = 50 ; % reference speed variation of the mixing matrix

%% Mixing matrix

typmel = 'nonstat';

switch typmel
    case 'stat'
        A = [1 0.75 -0.1; -0.5 1 0.3;-0.5 1 1];
        B = inv(A);
        z = A*y;
        for t=1:T
            heapA(:,:,t) = A;
            heapB(:,:,t) = B;
        end

    case 'nonstat'
        heapA = zeros(N,N,T);
        heapB = zeros(N,N,T);
        z = zeros(N,T);
        for t=1:T
            A = [1+0.3*cos(5*pi*f0*t/T) 0.75+0.4*cos(3*pi*f0*t/T);
                -0.7+0.5*cos(2*pi*f0*t/T) 1+0.1*cos(8*pi*f0*t/T)];
            B = inv(A);
            z(:,t) = A*y(:,t);
            c(t)=cond(A);
            heapA(:,:,t) = A;
            heapB(:,:,t) = B;
        end
end


%% Analysis
wav_typ = 'sharp';
wav_param = 1000;

NbScales = 150;
scales = 2.^(linspace(1,6,NbScales));


t = linspace(0,(T-1)/Fs,T);
nu0 = Fs/4;
freqdisp = [8 4 2 1 0.5]; % Displayed frequencies in kHz
sdisp = log2(nu0./(1e3*freqdisp)); % corresponding log-scales

figure; title('Sources and observations');
for n=1:N
    Wy = cwt_JEFAS(y(n,:),scales,wav_typ,wav_param);
    subplot(2,N,n); imagesc(t,log2(scales),abs(Wy)); yticks(sdisp); yticklabels(freqdisp);
    ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
    Wz = cwt_JEFAS(z(n,:),scales,wav_typ,wav_param);
    subplot(2,N,N+n); imagesc(t,log2(scales),abs(Wz)); yticks(sdisp); yticklabels(freqdisp);
    xlabel('Time (s)'); ylabel('Frequency (kHz)'); colormap(flipud(gray)); set(gca,'fontsize',18);
end