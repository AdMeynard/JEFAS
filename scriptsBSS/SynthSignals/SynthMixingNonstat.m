addpath('../../cwt');
clear all; close all;

K = 20;
for k= 1:K
%% Stationary sources 
    N = 2;
    T = 44100;
    x = randn(T,N);

    [x1,Sx1] = BandPassApprox(x(:,1),floor(0.06*T),floor(0.1*T));
    [x2,Sx2] = BandPassApprox(x(:,1),floor(0.15*T),floor(0.2*T));

    Sx = [Sx1'; Sx2'];

    %% Time warpings

    [y1,~,dgamma1] = chirpwarp(x1,1);
    [y2,~,dgamma2] = sinewarp(x2,2,0.02);

    y = [y1';y2'];
    dgamma = [dgamma1'; dgamma2'];

    %% Mixing matrix

    typmel = 'nonstat';

    switch typmel
        case 'stat'
            A = [1 0.75; -0.5 1];
            z = A*y;

        case 'nonstat'
            heapA = zeros(N,N,T);
            heapB = zeros(N,N,T);
            z = zeros(N,T);
            for t=1:T
                A = [1+0.3*cos(5*pi*t/T) 0.75+0.4*cos(3*pi*t/T);
                    -0.5+0.5*cos(11*pi*t/T) 1+0.1*cos(8*pi*t/T)];
                B = inv(A);
                z(:,t) = A*y(:,t);
                c(t)=cond(A);
                heapA(:,:,t) = A;
                heapB(:,:,t) = B;
            end
            figure;
            plot(c);
            title('A Condition number');
            figure;
            b11=heapB(1,1,:); b11 = b11(:);
            b12=heapB(1,2,:); b12 = b12(:);
            b21=heapB(2,1,:); b21 = b21(:);
            b22=heapB(2,2,:); b22 = b22(:);
            plot(1:T,b11,'b',1:T,b12,'k',1:T,b21,'r',1:T,b22,'g');
            title('B entries');
        yK(:,:,k) = y;
        zK(:,:,k) = z;
        heapAK(:,:,:,k) = heapA;
        heapBK (:,:,:,k) = heapB;
        SxK (:,:,k) = Sx;
        dgammaK (:,:,k) = dgamma; 
    end
   
end

