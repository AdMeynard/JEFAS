clear all; close all; clc;

load ../scriptsIEEE_TASLP/results/res_wind.mat
addpath('../cwt/') ;

%% Synthesize a broadband locally time-warped signal 
Fs = 2^11 ;
T = Fs;
t = linspace(0,(T-1)/Fs,T);


padL = 30 ;
im = 181;
Sx = [zeros(1,padL) Sx(1:im)];
nS = length(Sx) ;
fS = (0:(nS-1))/(2*(nS-1)) ;
f = (0:(T/2))/T ;
Sx = interp1(fS,Sx,f) ;
Sx = [Sx fliplr( Sx(2:(end-1)) ) ] ;
x = randn(T,1) ;
X = fft(x) .* sqrt(Sx(:)) ;
x = ifft(X) ;
 

dgammaML = smooth(dgammaML,3000) ;
dgammaML = dgammaML(4500:end);

TG = length(dgammaML) ;
tG = linspace(0,(TG-1)/TG,TG);
dgamma = interp1(tG,dgammaML,t) ;
dgamma = dgamma/mean(dgamma) ;

gamma = cumsum(dgamma)/T ;

y0 = sigwarp( x,gamma,dgamma ) ;
 
sigmay = 2e-1 ;
y = y0(:) + sigmay*randn(T,1);

save('../signals/wind4JEFASS','Fs','y0','dgamma','sigmay','y','Sx') ;

%% Show CWT
t = linspace(0,(T-1)/Fs,T) ;

Ms = 1000 ;
fmin = floor(0.058*T) ;
fmax = floor(0.42*T) ;

wav_typ = 'sharp'; % wavelet type (cf. cwt_JEFAS.m)
wav_param = 100 ;
freqtick = [100 200 400 800];

figure;

subplot(5,2,1); 
plot(t,x,'k','linewidth',2) ;
ylabel('Signal');
xticklabels([]);
set(gca,'fontsize',18);

subplot(5,2,[3 9]); 
[Wx,~,~] = display_cwt_JEFAS(x,Fs,fmin,fmax,Ms,wav_typ,wav_param,freqtick) ;
colormap(1-gray) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',18);

subplot(5,2,2); 
plot(t,y,'k','linewidth',2) ;
ylabel('Signal');
xticklabels([]);
set(gca,'fontsize',18);

subplot(5,2,[4 10]); 
[Wy,t,s] = display_cwt_JEFAS(y,Fs,fmin,fmax,Ms,wav_typ,wav_param,freqtick) ;
colormap(1-gray) ;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'fontsize',18);

figure;
plot(t,dgamma,'k','linewidth',2);
xlabel('Time (s)'); ylabel("\gamma'");
set(gca,'fontsize',18);