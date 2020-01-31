function yt = synthesis(W, MatPsi)
% Synthesize the signal from its adapted transform

MatfPsi = fft(MatPsi,[],1);

hatW = fft(W,[],2);
yt = sum(ifft(hatW.*(MatfPsi.'),[],2));
yt = real(yt);