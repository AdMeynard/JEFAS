function [SxSp,freqMax] = sparsifySpectrum(Sx,thres)

TS = length(Sx) ;
freqnorm = (0:(TS-1))/TS;

tmp = islocalmax(Sx) ;

SxSp = zeros(length(Sx),1) ;
SxSp(tmp) = Sx(tmp) ;
SxSp(SxSp<thres*max(SxSp)) = 0 ;

tmp(SxSp<thres*max(SxSp)) = 0 ;
tmp(freqnorm>.5) = 0 ; % remove negative frequencies
freqMax = freqnorm(tmp) ;