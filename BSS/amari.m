function indexdB = amari(Best, A)
% Amari index

G = Best*A; % should be close to identity
N = size(G,1);

index = (1/( 2*N*(N-1) ))*( sum( sum(abs(G).^2,2)./max(abs(G).^2,[],2) ) + ...
    sum( sum(abs(G).^2,1)./max(abs(G).^2,[],1) ) - 2*N );

indexdB = 10*log10(index);