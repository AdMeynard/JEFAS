function [y,gamma,dgamma] = newarp(x,f_w,Index)
% SINEWARP:   warp a vector using sinusoidal warping function derivative
%
% usage:    [y,gamma,dgamma] = sinewarp(x,f1,Index)
%
% Input:
%   x: original signal
%   f1: frequency of warping function
%   Index: amplitude of warping function
%
% Output:
%   y: output signal
%   gamma: warping function
%   dgamma: warping function derivative
%
% Conventions:
%   signals and warping functions are column vectors

if (Index <=0)
    error('Index is too large')
end

N = length(x);

ts = linspace(0,1,N)';

log2_dgamma = sin(2*pi*f_w*ts).*exp(-Index*ts); %sin(2*pi*f_w*ts).*exp(-Index*ts);%ex-cos %Index*atan((ts-0.5)*10);
dgamma = 2.^log2_dgamma;
mu = mean(dgamma);
dgamma = dgamma/mu;

gamma = zeros(N,1);
for k = 2:N
    gamma(k) = gamma(k-1) + dgamma(k);
end
gamma = gamma/gamma(N);

% y0 = sqrt(dgamma).*interp1(ts,x,gamma,'spline');
% y = y0';

y = interp1(ts,x,gamma,'spline');
end