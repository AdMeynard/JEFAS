function [pileA0, pileB0] = sobi_nonstat(z,vectau)
[N,T] = size(z);
dtau = vectau(2) -vectau(1);
ltau = length(vectau);

pileA0 = zeros(N,N,ltau);
pileB0 = zeros(N,N,ltau);
k = 1;
for tau=vectau
   Tau = tau:min((tau+dtau-1),T); 
   ztau = z(:,Tau);
   pileA0(:,:,k) = sobi(ztau,N);
%    pileB0(:,:,k) = inv(pileA0(:,:,k));
    B0 = inv(pileA0(:,:,k));
    pileB0(:,:,k) = B0./sqrt(sum(abs(B0).^2,2)); % rows normalization
   k = k+1;
end