function C = calc_cov_sparse(scales,S,theta)
% CALC_COV_SPARSE Computation the covariance matrix with sparse coefficents
% usage:	C = calc_cov_sparse(scales,S,theta)
%
% Input:
%   scales: vector of scales 
%   S: spectrum of X (row vector)
%   theta: Time warping parameter 
%
% Output:
%   C: covariance matrix

NbScales = length(scales);
T = length(S);
omega = (1:(T-1))*2*pi/T; % remove zero frequency for wavelets

scalesAct = 2^(-theta)*pi./(2*(omega(S(1:floor(T/2))~=0)+eps)); % scales where the signal is supposed to be active
k = 1;
for s = scalesAct
    [proxi,indices] = sort(abs(scales - s));
    ind(k,:) = indices(1:2);
    distances(k,:) = proxi(1:2);
    k = k+1;
end

d = zeros(NbScales,1);
d(ind(:,1)) = S(S(1:floor(T/2))~=0).*distances(:,2)./(distances(:,1)+distances(:,2));
d(ind(:,2)) = S(S(1:floor(T/2))~=0).*distances(:,1)./(distances(:,1)+distances(:,2));
C = diag(d);

end