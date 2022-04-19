function [SxSp,freqMax] = sparsifySpectrum(Sx,thres)
%SPARSIFYSPECTRUM	sparsify the power spectrum from its local maxima
% usage:	[SxSp,freqMax] = sparsifySpectrum(Sx,thres)
%
% Input:
%   Sx: non-sparse power spectrum
%   thres: sprasification threshold (0<thres<1)
% 
% Output:
%   SxSp : sparsified power spectrum
%   freqMax : active frequencies in the sparsified power spectrum


TS = length(Sx) ;
freqnorm = (0:(TS-1))/TS;

tmp = islocalmax(Sx) ;

SxSp = zeros(length(Sx),1) ;
SxSp(tmp) = Sx(tmp) ;
SxSp(SxSp<thres*max(SxSp)) = 0 ;

tmp(SxSp<thres*max(SxSp)) = 0 ;
tmp(freqnorm>.5) = 0 ; % remove negative frequencies
freqMax = freqnorm(tmp) ;