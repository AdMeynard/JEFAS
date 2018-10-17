function [Prefm, Prefw] = matrice_preference(haty,tau1,Deltat)
% tau1>= Deltat+1
% tau1<=T+1-Deltat
T = size(haty,2);

men = haty(:,max(tau1-Deltat,1):(tau1-1)).';
women = haty(:,tau1:min(tau1+Deltat-1,T)).';

R = matrice_corr(men, women,Deltat);
[~,Prefm] = sort(-R,2);
[~,Prefw] = sort(-R,1);
Prefw = Prefw.';

end
