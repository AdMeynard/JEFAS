function [W, MMSigmay] = SynthTSrep(y,sigmay,prior,scales,wav_typ,wav_param,Sx,theta,varargin) 
%SYNTHTSREP	generate the synthesis time-scale representation of a TW signal given the warping function and underlying spectrum
% usage:	[W, MMSigmay] = SynthTSrep(y,sigmay,prior,scales,wav_typ,wav_param,Sx,theta,TT,Delta) 
%
% Input:
%   y: signal
%   sigmay: noise variance
%   prior: shape of the a priori covariance matrix ('wavelet', 'sharp' or 'sparse')
%   scales:  vector of scales
%   wav_typ: wavelet type:
%     wav_typ='sharp': The sharp wavelet
%     wav_typ='dgauss': analytic derivative of gaussian
%   wav_par: parameter depending on the wavelet type
%     if wav_typ='sharp': wav_par = -ln(epsilon)>0 where epsilon=value of \hat{\psi}(Fs/2)
%     if wav_typ='dgauss': wav_par = number of vanishing moments
%   Sx: underlying spectrum
%   theta: WP parameters
%   TT: slicing size for W
%   Delta: overlap
% 
% Output:
%   W: adapted represnetation
%   MMSigmay: basis for signal covariance matrix estimation

T = length(y) ;

if isempty(varargin)
    TT = min(1024,round(T/2)) ;
else
    TT = varargin{1} ;
end

if length(varargin)~=2
    Delta = 0.75*TT ;
else
    Delta = varargin{2} ;
end


[M_psi, ~] = bas_calc_dcov(scales,wav_typ,wav_param,T);
MatPsi = ifft(M_psi.',[],1); % Synthesis coefficients

if strcmp(prior,'sharp')
    sigma_sharp = log2(scales(2)/scales(1)) ;
    [M_psi, ~] = bas_calc_dcov_sharp(scales,sigma_sharp) ;
end

T = length(y);
TS = length(Sx);
omega = (0:(TS-1))*2*pi/TS;

delta = floor( (TT-Delta)/2 );
K = ceil( (T-TT)/Delta ) + 1; % number of blocks

for dec=0:(TT-1)
    Mdec = circshift(MatPsi,dec,1);
    MN{dec+1} = Mdec(1:TT,:);
end

nMax = 0;
for k = 1:K % k-th block
    nn = ( (k-1)*Delta+1 ):( (k-1)*Delta+TT );
    nn = nn((nn>=1)&(nn<=T)); % time of the k-th block
    NN = length(nn);
    ycourt = y(nn); % we divide the signal in segments of size TT
    
    if isreal(MatPsi)
        Sigmay = sigmay^2*eye(NN);
    else
        Sigmay = 2*sigmay^2*eye(NN);
    end
    
    for n = nn
        Mn = MN{n-min(nn)+1};
        if NN ~= TT
            Mn = Mn(1:NN,:);
        end
        if n > nMax % avoid redundancies
            switch prior 
                case 'wavelet'
                    C{n} = calc_cov_synth(M_psi,Sx,omega,TS,theta(n));
                case 'sharp'
                    C{n} = calc_cov_synth_sharp(M_psi,Sx,omega,scales,theta(n));
                case 'sparse'
                    C{n} = calc_cov_sparse(scales,Sx,theta(n));
            end
        end
        CMn{n-min(nn)+1} = C{n}*Mn';
        Sigmay = Sigmay + Mn*CMn{n-min(nn)+1};
    end
    Sigmay = real(Sigmay);
    MMSigmay{k} = Sigmay; % we store successive Sigmay matrices
    
    % Adapted representation computation
    invSy = Sigmay\ycourt;
    if k == 1 % first block
        for n = nn
            W(:,n) = CMn{n-min(nn)+1} * invSy ;
        end
    elseif k<K
        nnC = nn( delta:(TT - delta - 1) );
        for n = nnC
            W(:,n) = CMn{n-min(nn)+1} * invSy ;
        end
    else
        nnC = nn( delta:end );
        for n = nnC
            W(:,n) = CMn{n-min(nn)+1} * invSy ;
        end
    end
    
    nMax = nn(end);
end
