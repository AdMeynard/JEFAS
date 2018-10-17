function [N, M, ptszd, ptsd] = selecpts(z,pp,eps1,eps2,eps3,eps4)

[N,T] = size(z); % original sig
tt = 1:pp:T; % subsampling in order to have correct size matrices
TT = length(tt); % subsampledsig

% Construction des RTFQ
Dz = zeros(TT,TT,N,N);
for n = 1:N
    for m = 1:N
        Dznm = tfrpwv([hilbert(z(n,:)).', hilbert(z(m,:)).'],tt,TT);
        Dz(:,:,n,m) = Dznm;
    end
end
Mat = permute(Dz,[3 4 1 2]);
Mat = reshape(Mat,N,N,TT^2);

% Sélection des points
iDz = sum(sum(abs(imag(Dz)).^2,4),3);
rDz = sum(sum(abs(real(Dz)).^2,4),3);
% ZD
ptszd = find((iDz > eps1)&(iDz > eps2*rDz));
N = Mat(:,:,ptszd);
% D
ptsd = find((iDz < eps3)&(rDz > eps4));
M = Mat(:,:,ptsd);
