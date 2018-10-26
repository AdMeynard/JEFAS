function [y,gamma,dgamma] = sinewarp(x,f1,Index)
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

if (Index >=1/(2*pi*f1))
    error('Index is too large')
end

N = length(x);

ts = linspace(0,1,N)';

gamma = ts + Index*sin(2*pi*f1*ts);

dgamma = 1 + 2*pi*f1*Index*cos(2*pi*f1*ts);

% y0 = sqrt(dgamma).*interp1(ts,x,gamma,'spline');
% y = y0';

y = sigwarp( x,gamma,dgamma );
end

