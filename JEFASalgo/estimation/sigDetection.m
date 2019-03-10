function act = sigDetection(y,ratio,scales,wav_typ,wav_param)

T = length(y);
act = ones(1,T);

Wy = cwt_JEFAS(y,scales,wav_typ,wav_param);

evolPower = sum( abs(Wy).^2 );
medPower = median(evolPower);
threshold = ratio*medPower;

act(evolPower <= threshold) = 0;