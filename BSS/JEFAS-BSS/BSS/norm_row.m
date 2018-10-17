function [c, ceq] = norm_row(B)
P = size(B,1);

v = sqrt(sum(abs(B).^2,2));
N = length(v);
ceq = v - ones(P,1);
c=[];