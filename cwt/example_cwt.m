clear all; close all; clc;
addpath('../signals');

load('wolf_howling');
T = length(y);
t = linspace(0,(T-1)/Fs,T);

scales = 2.^linspace(3,6,200);
wav_typ = 'sharp';

wav_par1 = 50;
W1 = cwt(y,scales,wav_typ,wav_par1);
wav_par2 = 500;
W2 = cwt(y,scales,wav_typ,wav_par2);

figure;
subplot(1,2,1);imagesc(t,scales,abs(W1)); 
set(gca,'yticklabel',[]); xlabel('Time (s)');
subplot(1,2,2);imagesc(t,log2(scales),abs(W2));

nu0 = Fs/4;
sobs = cellfun(@str2num,get(gca,'yticklabel'));
fobs = round(nu0./2.^sobs);
set(gca,'yticklabel',fobs);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); colormap(flipud(gray));


nu1 = Fs/2;

N = 100001;
fsw = 1e6;
nu = linspace(0,fsw,N);

hatpsi1 = exp( -wav_par1*( 0.5*(nu/nu0 + nu0./nu) - 1 ) );
hatpsi2 = exp( -wav_par2*( 0.5*(nu/nu0 + nu0./nu) - 1 ) );
figure;
plot(nu,hatpsi1,'r',nu,hatpsi2,'b','linewidth',2); grid on; axis tight; xlim([0 Fs/2]);
legend('\epsilon = 10^{-50}','\epsilon = 10^{-500}');
xlabel('Frequency (Hz)'); ylabel({'$\hat\psi_\sharp$'},'interpreter','latex');

t = linspace(-(N-1)/(2*fsw),(N-1)/(2*fsw),N);
psi1 = ifftshift(ifft(hatpsi1));
psi2 = ifftshift(ifft(hatpsi2));
figure;
subplot(2,1,2);plot(t,real(psi1),'r',t,imag(psi1),'r:','linewidth',2); xlim([-0.001 0.001]); grid on;
legend('Real part','Imaginary part');
xlabel('Time (s)'); ylabel({'$\psi_\sharp$'},'interpreter','latex');
subplot(2,1,1); plot(t,real(psi2),'b',t,imag(psi2),'b:','linewidth',2); xlim([-0.001 0.001]); grid on;
legend('Real part','Imaginary part');
ylabel({'$\psi_\sharp$'},'interpreter','latex');