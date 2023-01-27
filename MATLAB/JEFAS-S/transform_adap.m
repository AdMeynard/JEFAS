function [W, MMSigmay] = transform_adap(y,sigmay,priorList,TT,Delta,M_psi,Sx,theta,MatPsi) 
%TRANSFORM_ADAP	determine the adapted transform of a signal given the warping function and underlying spectrum
% usage:	[W, MMSigmay] = transform_adap(y,sigmay,priorList,TT,Delta,M_psi,Sxest,theta,MatPsi) 
%
% Input:
%   y: signal
%   sigmay: noise variance
%   priorList: cell of shape {prior, scales} where
%     prior: shape of the a priori covariance matrix ('wavelet', 'sharp' or 'sparse')
%     scales: vector of scales (only when prior~='wavelet' because scales is included in M_psi)
%   TT: slicing size for W
%   Delta: overlap
%   M_psi: base for covariance computation
%   Sx: underlying spectrum
%   theta: WP parameters
%   MatPsi: characterize the transform shape (wavelet)
% 
% Output:
%   W : adapted transform
%   MMSigmay : basis for signal covariance matrix estimation

prior = priorList{1} ;
if ~strcmp(prior,'wavelet')
    scales = priorList{2} ;
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