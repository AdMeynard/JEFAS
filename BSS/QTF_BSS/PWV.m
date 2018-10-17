function pvwxy = PWV(sig1,sig2,h)

siglength = length(sig1);

if isreal(sig1)
    sig1 = hilbert(sig1(:));
else
    sig1 = sig1(:);
end

if isreal(sig1)
    sig2 = hilbert(sig2(:));
else
    sig2 = sig2(:);
end

v1 = [sig1(end) zeros(1,siglength-1)]; 
% v1 = circshift(sig1,1); % périodique
M1 = hankel(sig1,v1);

v2 = [sig2(1)'; zeros(siglength-1,1)];
% v2 = [sig2(1)'; flipud(conj(sig1(2:end)))]; % périodique
M2 = toeplitz(v2,sig2');

Mh = repmat(h(:),1,siglength);

fpvw = 2*M1.*M2.*Mh;
pvwxy = fft(fpvw,[],1);
    
end