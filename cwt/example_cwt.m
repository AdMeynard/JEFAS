clear all; close all; clc;
addpath('../signals');

load('wolf_howling');
T = length(y);
t = linspace(0,(T-1)/Fs,T);

%% Two wavelet transforms of the same signal

scales = 2.^linspace(3,6,200);
wav_typ = 'sharp';

wav_par1 = 50;
W1 = cwt(y,scales,wav_typ,wav_par1); % wavelet transform
wav_par2 = 500;
W2 = cwt(y,scales,wav_typ,wav_par2);

nu0 = Fs/4; % wavelet central frequency
freqdisp = [1.50 1.25 1.00 0.75 0.50 0.25]; % Displayed frequencies
sdisp = log2(nu0./(1e3*freqdisp)); % coreesponding log-scales

figure; colormap(flipud(gray))
subplot(1,2,1);imagesc(t,log2(scales),log1p(abs(W1)/0.1)); % scalogram
xlabel('Time (s)'); ylabel('Frequency (kHz)');
yticks(sdisp); yticklabels(freqdisp);
set(gca,'fontsize',24);

subplot(1,2,2);imagesc(t,log2(scales),log1p(abs(W2)/0.1));
xlabel('Time (s)'); ylabel('Frequency (kHz)');
yticks(sdisp); yticklabels(freqdisp);
set(gca,'fontsize',24);

yrec1 = icwt(W1,scales,wav_typ,wav_par1); % reconstructed signal from W1
yrec1 = std(y)*yrec1/std(yrec1); % normalize the reconstructed signal (because ICWT is given up to a constant)
yrec2 = icwt(W2,scales,wav_typ,wav_par2);
yrec2 = std(y)*yrec2/std(yrec2);

figure; 
subplot(2,1,1);plot(t,y,'r',t,yrec1,'b',t,yrec2,'c--','linewidth',2);
title('Wavelet transform inversion'); xlabel('Time (s)'); ylabel('Signals'); legend('Original signal', 'Signal reconstructed from W1','Signal reconstructed from W2');
subplot(2,1,2);plot(t,y,'r',t,yrec1,'b',t,yrec2,'c--','linewidth',2); xlim([0.35, 0.45]); % zoom
title('Zoom'); xlabel('Time (s)'); ylabel('Signals'); legend('Original signal', 'Signal reconstructed from W1','Signal reconstructed from W2');

%% Plot the wavelet and its Fourier transform

nu1 = Fs/2; % wavelet frequency such that hatpsi(nu1) = epsilon

N = 100001;
fsw = 1e6;
nu = linspace(0,fsw,N);

hatpsi1 = exp( -wav_par1*( 0.5*(nu/nu0 + nu0./nu) - 1 )/(0.5*(nu1/nu0 + nu0./nu1) - 1 ) ); % Fourier transform of psi
hatpsi2 = exp( -wav_par2*( 0.5*(nu/nu0 + nu0./nu) - 1 )/(0.5*(nu1/nu0 + nu0./nu1) - 1 ) );
figure;
plot(nu,hatpsi1,'r',nu,hatpsi2,'b','linewidth',2); grid on; axis tight; xlim([0 Fs/2]);
legend('\epsilon = e^{-50}','\epsilon = e^{-500}');
xlabel('Frequency (Hz)'); ylabel({'$\hat\psi_\sharp$'},'interpreter','latex');

t = linspace(-(N-1)/(2*fsw),(N-1)/(2*fsw),N);
psi1 = ifftshift(ifft(hatpsi1)); % the sharp wavelet !
psi2 = ifftshift(ifft(hatpsi2));
figure;
subplot(2,1,2);plot(t,real(psi1),'r',t,imag(psi1),'r:','linewidth',2); xlim([-0.002 0.002]); grid on;
legend('Real part','Imaginary part');
xlabel('Time (s)'); ylabel({'$\psi_\sharp$'},'interpreter','latex');
subplot(2,1,1); plot(t,real(psi2),'b',t,imag(psi2),'b:','linewidth',2); xlim([-0.002 0.002]); grid on;
legend('Real part','Imaginary part');
ylabel({'$\psi_\sharp$'},'interpreter','latex');