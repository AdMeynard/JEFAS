function R = matrice_corr(men, women,Dtau)
Tm = size(men,1);
Tw = size(women,1);

Lm = linspace(0,1,Tm);
Lw = linspace(0,1,Tw);
L = linspace(0,1,Dtau);


tf_men = abs(fft(men,[],1)).^2; % spectrum
tf_men = tf_men./sqrt(sum(abs(tf_men).^2)); % normalization
tf_men = interp1(Lm,tf_men,L); % vector length = Dtau

tf_women = abs(fft(women,[],1)).^2;
tf_women = tf_women./sqrt(sum(abs(tf_women).^2));
tf_women = interp1(Lw,tf_women,L);

% A voir
% R = zeros(N);
% for n=1:N
%     tf_man = tf_men(:,n);
%     for k=1:N
%         tf_woman = tf_women(:,k);
%         R(n,k) = max(abs(xcorr(tf_man,tf_woman)));
%     end
% end

R = tf_men'*tf_women; % correlation matrix