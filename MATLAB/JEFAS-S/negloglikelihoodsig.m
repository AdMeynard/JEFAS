function nllk = negloglikelihoodsig(y,Sigmay)
%NEGLOGLIKELIHOODSIG negative log likelihood of the signal
% usage:	nllk = negloglikelihoodsig(y,Sigmay)
%
% Input:
%   y: signal
%   Sigmay: signal covariance matrix
%
% Output:
%   nllk: value of the negative log likelihood

nllk = y.' * (Sigmay\y) + real( logdet(Sigmay) ); % Should decrease through iterations in EM algorithm