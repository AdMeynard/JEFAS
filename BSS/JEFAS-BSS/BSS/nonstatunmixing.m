function haty = nonstatunmixing(z,pileB,vectau)

Sizetau = length(vectau);
T = size(z,2);

% piecewise umixing
for k = 1:(Sizetau-1)
    tau0 = vectau(k);
    tau1 = vectau(k+1) - 1;
    haty(:,tau0:tau1) = pileB(:,:,k)*z(:,tau0:tau1); % k-th segment
end

tau0 = vectau(Sizetau);
haty(:,tau0:T) = pileB(:,:,Sizetau)*z(:,tau0:T); % last segment

