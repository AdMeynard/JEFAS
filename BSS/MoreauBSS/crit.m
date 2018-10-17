function c = crit(U,N,M)
% Critere a minimiser

K = size(M,3);
P = size(N,3);
m = min(K,P);

c = 0;
for k = 1:m
    c = c - sum(abs(diag(U'*M(:,:,k)*U)).^2) + sum(abs(diag(U'*N(:,:,k)*U)).^2);
end

if K>P
    for k = (m+1):K
       c = c - sum(abs(diag(U'*M(:,:,k)*U)).^2);
    end
elseif K<P
    for k = (m+1):P
       c = c + sum(abs(diag(U'*N(:,:,k)*U)).^2);
    end
end


