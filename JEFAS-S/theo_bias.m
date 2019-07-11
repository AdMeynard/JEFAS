function bias = theo_bias(y0,TT,Delta,MMSigmay,sigmay)
% Bias on the sythesis of the signal from the adapted transform

T = length(y0);
delta = floor( (TT-Delta)/2 );
K = ceil( (T-TT)/Delta ) + 1; % nb blocs
for k = 1:K
    nn = ( (k-1)*Delta+1 ):( (k-1)*Delta+TT );
    nn = nn((nn>=1)&(nn<=T)); % instants du bloc k
    NN = length(nn);
    y0k = y0(nn);
    b = -2*sigmay^2*(MMSigmay{k}\y0k);

    interv = delta:(TT - delta - 1);
    if k == 1 % premier bloc
        bias(nn) = b;
    elseif k<K
        bias( nn(interv) ) = b(interv);
    else
        bias( nn(delta:end) ) = b(delta:end);
    end
end

