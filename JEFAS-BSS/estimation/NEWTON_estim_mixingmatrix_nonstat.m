function heapBoptim = NEWTON_estim_mixingmatrix_nonstat(A0,Wz,vectau,M_psi,Sx,dgamma,L)

[Ms,~,N] = size(Wz);
Kmat = length(vectau);

heapBoptim = zeros(N,N,Kmat);
a = A0(:);
A = A0;

% errmat = eps_bss*ones(N); % prevent brutal variations of the unmixing matrix

delta = 1:N;
k = 1;
for tau0 = vectau
    wz = reshape(Wz(:,tau0,:),Ms,N).';
    % covariance matrices:
    for n=1:N
        C(:,:,n) = calc_cov(M_psi,Sx(n,:),log2(dgamma(n,tau0)));
        M(:,:,n) = real(wz * (C(:,:,n) \ wz') ) ;
    end
    
    % log-likelihood optimization:
    for l=1:L
        aold = a;
        Aold = A;
        
        J = zeros(N^2, N^2);
        for n=1:N
            Mna = M(:,:,n) \ Aold ;
            ga( ((n-1)*N+1) : n*N ) = Aold(:,n).' * Mna - double(delta==n) ;
            J( ((n-1)*N+1) : n*N, ((n-1)*N+1) : n*N ) = Mna.' ;
            J( ((n-1)*N+1) : n*N, : )  = J( ((n-1)*N+1) : n*N, : ) + kron( eye(N), Mna(:,n).' );
        end  
        a = aold - J\ga(:) ;
        A = reshape(a,N,N);
    end
    heapBoptim(:,:,k) = inv(A);
    
    k = k+1;
end