function y = sigwarp( x,gamma,dgamma )
%WARP unitary signal warping
%   usage:  y = warp( x,gamma,dgamma )

N = length(x);
ts = linspace(0,1,N);
y = sqrt(dgamma).*interp1(ts,x,gamma,'spline');
end

