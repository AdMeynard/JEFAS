function [y,Sx] = BandPassApprox(x,N_l,N_h)
% LowPassApprox -- Low pass approximation
%  Usage
%    [y,err] = LowPassApprox(x,N)
%  Input Variables:
%    x: input signal(column vector)
%    N: cutoff frequency
%  Output Variables:
%    y:  low pass approximation
%    err: approximation error
%

L = length(x);

% Verification
if ~iscolumn(x)
	error('x must be a column vector')
end
if N_h>L/2
	error('N must be smaller than L/2');
end

% Calculation
xchap = fft(x);
mask = zeros(L,1);
vec = N_l:N_h;
if ~mod(L,2)
	mask(N_l:N_h)= sqrt(1 - cos(2*pi*(vec - N_l)/(N_h - N_l)));
	mask((L-N_h+2):(L-N_l+2))=mask(N_l:N_h);
else
	mask(N_l:N_h)= sqrt(1 - cos(2*pi*(vec - N_l)/(N_h - N_l)));
	mask((L-N_h):(L-N_l))=mask(N_l:N_h);
end
y = ifft(xchap.*mask);

y = real(y);

% err = norm(y-x)/norm(x);
Sx = mask.^2;

end