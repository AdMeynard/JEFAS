clear  all; close all; clc;
addpath('../cwt/') ;

Fs = 2^10 ;
T = Fs ;

%% Generate a rapidly varying warping function
s1 = load('../signals/SingleComp.mat') ;
SigLength0 = length(s1.s) ;
ind0 = linspace(0,1,SigLength0) ;

dgamma = s1.if1/mean(s1.if1) ;
ind1 = linspace(0,1,T) ;
dgamma = interp1(ind0,dgamma,ind1) ;
gamma = cumsum(dgamma)/Fs ;

%% Synthesize a broadband locally time-warped signal 
t = linspace(0,(T-1)/Fs,T) ;

f1  = round(0.18*Fs);
phi1 = 2*pi*rand ;
a1 = 1 ;
x1 = a1*cos(2*pi*f1*t+ phi1) ;
f2  = round(0.15*Fs);
phi2 = 2*pi*rand ;
a2 = 2 ;
x2 = a2*cos(2*pi*f2*t+phi2) ;
x = x1 + x2 ;

Sx = zeros(T,1);
Sx(f1+1) = a1^2;
Sx(T-f1) = a1^2;
Sx(f2+1) = a2^2;
Sx(T-f2) = a2^2;

sigmay = 5e-2;
y0 = sigwarp( x,gamma,dgamma )' ;
y = y0 + sigmay*randn(T,1);

save('../signals/TwoSineWaves_fast','Fs','Sx','y0','dgamma','sigmay','y');

%% Show CWT

Ms = 1000 ;
fmin = floor(0.048*T) ;
fmax = floor(0.40*T) ;

wav_typ = 'sharp'; % wavelet type (cf. cwt_JEFAS.m)
wav_param = 500 ;
freqtick = [50 100 200 300 400];

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