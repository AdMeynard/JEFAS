function [W, MMSigmay] = transform_adap_sparse(y,sigmay,TT,Delta,scales,Sxest,thetaE,MatPsi) 
% Calcul de W

T = length(y);
Sxest = Sxest(:);


delta = floor( (TT-Delta)/2 );
K = ceil( (T-TT)/Delta ) + 1; % nb blocs

for dec=0:(TT-1)
    Mdec = circshift(MatPsi,dec,1);
    MN{dec+1} = Mdec(1:TT,:);
end

nMax = 0;
for k = 1:K % k-eme bloc temporel
    nn = ( (k-1)*Delta+1 ):( (k-1)*Delta+TT );
    nn = nn((nn>=1)&(nn<=T)); % instants du bloc k
    NN = length(nn);
    ycourt = y(nn); % on saucissonne le signal en blocs de taille TT
    
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
        if n > nMax % evite les redondances
            C{n} = calc_cov_sparse(scales,Sxest,thetaE(n));
        end
        CMn{k,n} = C{n}*Mn';
        Sigmay = Sigmay + Mn*CMn{k,n};
    end
    Sigmay = real(Sigmay);
    MMSigmay{k} = Sigmay; % on stocke les matrice Sigmay succesives
    
    % Calcul de W
    invSy = Sigmay\ycourt;
    if k == 1 % premier bloc
        for n = nn
            W(:,n) = CMn{k,n} * invSy ;
        end
    elseif k<K
        nnC = nn( delta:(TT - delta - 1) );
        for n = nnC
            W(:,n) = CMn{k,n} * invSy ;
        end
    else
        nnC = nn( delta:end );
        for n = nnC
            W(:,n) = CMn{k,n} * invSy ;
        end
    end
    
    nMax = nn(end);
end