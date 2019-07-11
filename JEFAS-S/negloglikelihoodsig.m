function nllk = negloglikelihoodsig(y,Sigmay)
% Should decrease through iterations in EM algorithm

nllk = y.' * (Sigmay\y) + real( logdet(Sigmay) );