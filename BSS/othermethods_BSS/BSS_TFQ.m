function [A] = BSS_TFQ(z,pp, eps3, eps4, nn)

N = size(z,1);
[~, Dmat] = selecpts(z, pp, 0, 1, eps3, eps4);

Kd = size(Dmat,3);
for k= 1:Kd
    Dx = real(Dmat(:,:,k));
    [U,~,~] = svd(Dx);
    colA(:,k) = U(:,1);
end

[idx,C] = kmeans(colA',nn); % clustering
nC = histc(idx(:),1:nn); % frequency of each class
for k=1:N
    A(:,k) = C(nC==max(nC),:).';
    nC(nC==max(nC)) = 0;
end