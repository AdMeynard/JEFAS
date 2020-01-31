

%% Synthesize a broadband locally time-warped signal 
Fs = 2^13; %Fs = 2^10;
T = Fs;
t = linspace(0,(T-1)/Fs,T);

x = randn(T,1);
[x1,Sx1] = BandPassApprox(x,floor(0.10*T),floor(0.12*T));
[x2,Sx2] = BandPassApprox(x,floor(0.05*T),floor(0.08*T));
x = x1+x2;
Sx = Sx1 + Sx2;

% [y0,~,dgamma] = sinewarp(x,2,0.02);
[y0,gamma,dgamma] = newarp(x,2,2);

sigmay = 5e-2;
y = y0 + sigmay*randn(T,1);

save('../signals/toysigJEFASS','Fs','Sx','y0','dgamma');
save('../signals/toysig_short','Fs','Sx','y0','dgamma','sigmay','y');