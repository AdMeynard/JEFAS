function ordered_newsig = reordersig(oldsig,newsig)
N = size(newsig,1);

for n=1:N
    x = oldsig(n,:);
    for k=1:N
        y = newsig(k,:);
        M(n,k) = max(abs(xcorr(x,y)));
    end
end

[~,Prefm] = sort(-M,2);
[~,Prefw] = sort(-M,1);
Prefw = Prefw.';

stablematch = galeshapley(N,Prefm,Prefw);
ordered_newsig = newsig(stablematch,:);
