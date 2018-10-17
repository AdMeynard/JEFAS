function [A] = BSS_Moreau(Dmat,N, nn)

Kd = size(Dmat,3);
for k= 1:Kd
    Dx = real(Dmat(:,:,k));
    [U,~,~] = svd(Dx);
    colA(:,k) = U(:,1);
end

[idx,C] = kmeans(colA',nn); % clustering
nC = histc(idx(:),1:nn); % occurence de chaque classe
for k=1:N
    A(:,k) = C(nC==max(nC),:).'; % On ne prend que les 
    nC(nC==max(nC)) = 0;
end