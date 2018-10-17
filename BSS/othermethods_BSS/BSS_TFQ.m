function [A] = BSS_TFQ(Dmat,N, nn)

Kd = size(Dmat,3);
for k= 1:Kd
    Dx = real(Dmat(:,:,k));
    [U,~,~] = svd(Dx);
    colA(:,k) = U(:,1);
end

% Z = linkage(colA','ward');
% idx = cluster(Z,'Maxclust',6);

[idx,C] = kmeans(colA',nn); % clustering
nC = histc(idx(:),1:nn); % occurence de chaque classe
for k=1:N
    A(:,k) = C(nC==max(nC),:).';
%     u = (idx==find(nC==max(nC)));
%     A(:,k) = mean(colA(:,u),2);
    nC(nC==max(nC)) = 0;
end