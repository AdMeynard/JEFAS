function Sigmay = buildSigmay(MMSigmay,T,TT,Delta)
% Concatenate Small Sigmay to create full Sigmay


K = length(MMSigmay);

for k = 1:K % k-eme bloc temporel
    nn = ( (k-1)*Delta+1 ):( (k-1)*Delta+TT );
    nn = nn((nn>=1)&(nn<=T)); % instants du bloc k
    NN = length(nn);
    
    Sigmay(nn,nn) = MMSigmay{k} ;
end

% if k == 1 % premier bloc
%     Cy(nn,nn) = Sigmay ;
% elseif k<K
%     nnC = nn( delta:(TT - delta - 1) );
%     Cy(nnC,nnC) = Sigmay(nnC-(k-1)*Delta-1,nnC-(k-1)*Delta-1) ;
% else
%     nnC = nn( delta:end );
%     Cy(nnC,nnC) = Sigmay(nnC-(k-1)*Delta-1,nnC-(k-1)*Delta-1) ;
% end