function [WT, WSST, WSST2, fs, as, omega, omega2, tau, phipp, norm2psi] = Wsst2(s,gamma,mywav,nv)
% function Wsst : computes the synchrosqueezing transform of signal s.
% New formulae Cphi for reconstruction
%
% Inputs
%   s : input signal, power of 2
%   dt : sample period
%   gamma : threshold
%   mywav : string coding the mother wavelet (see mycwt)
%   nv : number of coefficient per octave
%   doplot : 0 or 1 for extern plotting
%
% Outputs : 
%   WT : the wavelet transform
%   WSST : the synchrosqueezed transform
%   WSST2: the second order
%   fs : frequency vector
%   as : scales vector


% nargin/out
if nargin<4
    nv = 32; % nbre coef    s par octave
    gamma = 0.01; % Seuil pour wavelet thresholding
    mywav = 'cmor2-1';
end

% Computing scales
s = s(:);
n = length(s);
noct = log2(n)-1;
na = noct*nv+1;
as = (2^(-1/nv)) .^ (0:1:na-1);

%noct = log2(n)-1;
%na = noct*nv;
%as = (2^(1/nv) .^ (1:1:na));

% Padding signal (symmetric)
[N, sx, n1] = mypad(s);

%% Wavelet transform and reassignment operators
WT = zeros(na, N); % CWT
omega = zeros(na, N); % Vertical operator
tau = zeros(na, N); % Horizontal one
phipp = zeros(na, N); % FM
omega2 = zeros(na, N); % Second-order IF
Dtau =  zeros(na, N);

sx = sx(:).';
xh = fft(sx);
%xi = (0:N-1)/N;
xi = (0:N-1)/2; % to change if signal not a power of 2

% Filter definition
if strncmp(mywav,'gmor',4)
    [v1 v2] = regexp(mywav,'[0-9.]*-[0-9.]');
    beta = str2num(mywav(v1:v2-2));
    gam = str2num(mywav(v2:end));
    filt = @(a) gmor(beta,gam,a*xi);
    func = @(x) 2*(exp(1)*gam/beta)^(beta/gam) * x.^beta .* exp(-x.^gam);
elseif strncmp(mywav,'cmor',4)
    [v1 v2] = regexp(mywav,'[0-9.]*-[0-9.]');
    Fb = str2num(mywav(v1:v2-2));
    Fc= str2num(mywav(v2:end));
    filt = @(a) cmor(Fb,Fc,a*xi);
    func = @(x) sqrt(Fb) * exp(-Fb^2*pi*(x-Fc).^2);
elseif strncmp(mywav,'bump',4)
    [v1 v2] = regexp(mywav,'[0-9.]*-[0-9.]');
    mu = str2num(mywav(v2:end));
    sigma = str2num(mywav(v1:v2-2));
    filt = @(a) bump(mu,sigma,a*xi);
    func = @(x) exp(1-1./(1-((x-mu)./sigma).^2));
end

% for each octave
norm2psi = zeros(1,length(as));
for ai = 1:length(as)
    a = as(ai);
    [psih, dfilt,~] = filt(a);
    norm2psi(ai) = norm(psih);
    % temporary values
    Wtmp = ifft(conj(psih).* xh); 
    Wnu = ifft(conj(a*xi.*psih) .* xh);
    Wnunu = ifft(conj((a*xi).^2.*psih) .* xh);
    Wp = ifft(conj(dfilt).*xh);
    Wnup = ifft(conj(a*xi.*dfilt).*xh);
    
    % store results
    WT(ai,:) = Wtmp;
    omega(ai,:) = 1/a*Wnu./Wtmp;
    tau(ai, :) = a*Wp./Wtmp/2/1i/pi;
    phipp(ai,:) = 2*1i*pi/a^2*(Wnunu.*Wtmp-Wnu.^2)./(Wtmp.^2+Wnup.*Wtmp-Wp.*Wnu);
    Dtau(ai,:) = (Wtmp.^2+Wnup.*Wtmp-Wp.*Wnu)./(Wtmp.^2);
end

%phipp(abs(real(tau)*n)<1)=0;
omega2 = real(omega - phipp.*tau);
omega = real(omega);
tau= real(tau);


WT = WT(:, n1+1:n1+n);
omega = omega(:, n1+1:n1+n);
tau = tau(:, n1+1:n1+n);
phipp = phipp(:, n1+1:n1+n);
omega2 = omega2(:, n1+1:n1+n);
Dtau = Dtau(:, n1+1:n1+n);

% figure(); 
% imagesc((0:n-1)/n,1./as,abs(Dtau)); set(gca,'YDir','normal');
% 
% pause
 
%%
omega(omega<0) = NaN;
omega(abs(WT)<gamma) = NaN;
tau(abs(WT)<gamma) = NaN;
omega2(omega2<0) = NaN;
omega2(abs(WT)<gamma) = NaN;
phipp(abs(WT)<gamma) = NaN;

%% Synchrosqueezing and reassignment
fs = 1./as;
WSST = zeros(size(WT));
WSST2 = zeros(size(WT));

for b=1:n
    for ai=1:length(as)
        if (~isnan(omega(ai,b)) && abs(WT(ai,b))> 2*gamma*sqrt(2)*norm2psi(ai)/sqrt(2*n)) %module 
            %if (abs(Dtau(ai,b)) > 0)
             k = round(1+nv*log2(omega2(ai,b)));
             if k>0 && k<=na 
              WSST2(k,b) = WSST2(k, b) + WT(ai, b); % WSST2
             end
%             else
%              k = round(1+nv*log2(omega(ai,b)));
%              if k>0 && k<=na
%             	WSST2(k,b) = WSST2(k, b) + WT(ai, b); % WSST
%              end
%             end 
            k = round(1+nv*log2(omega(ai,b)));
            if k>0 && k<=na
            	WSST(k,b) = WSST(k, b) + WT(ai, b); % WSST
            end
        end
    end % for ai...
end % for b

WSST = WSST/nv*log(2);
WSST2 = WSST2/nv*log(2);
end

