function [W, MMSigmay] = transform_adap(y,sigmay,TT,Delta,M_psi,Sxest,thetaE,MatPsi) 
% Calcul de W

T = length(y);
TS = length(Sxest);
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
    ycourt = y(nn); % we dice the signal in segments of size TT
    
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
            C{n} = calc_cov_synth(M_psi,Sxest,omega,TS,thetaE(n));
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
