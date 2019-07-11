clear all; close all;

N = 2^13;
t = linspace(0,(N-1)/N,N);

a = 200;
b1 = 500;
phi1 = a*t.^2/2 + b1*t;
b2 = 550;
phi2 = a*t.^2/2 + b2*t;
b3 = 650;
phi3 = a*t.^2/2 + b3*t;

x1 = cos(2*pi*phi1);
x2 = cos(2*pi*phi2);
x3 = cos(2*pi*phi3);

y12 = x1 + x2;
y13 = x1 + x3;
Ms = 200;
freqtick = [1000 700 500];

figure;
subplot(1,2,2);
Wy12 = display_cwt_JEFAS(y12,N,400,1100,Ms,'sharp',220,freqtick);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); colormap(1-gray);
set(gca,'FontSize',20);
subplot(1,2,1);
Wy13 = display_cwt_JEFAS(y13,N,400,1100,Ms,'sharp',220,freqtick);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); colormap(1-gray);
set(gca,'FontSize',20);