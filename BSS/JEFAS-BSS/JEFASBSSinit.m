function [pileB_init, vectau_init] = JEFASBSSinit(z, init_meth, varargin)

switch init_meth
    case 'psobi'
        Dt = varargin{1};
        vectauns = 1:Dt:T;
        [~, pile_Bsobins] = sobi_nonstat(z,vectauns);
        vectau_init = vectauns;
        pileB_init = pile_Bsobins;
    case 'sobi'
        N = size(z,1);
        [Asobi,~] = sobi(z,N,1000);
        B0 = inv(Asobi);
        B0 = B0./sqrt(sum(abs(B0).^2,2));
        vectau_init = 1;
        pileB_init = B0;
    case 'Moreau'
        N = size(z,1);
        eps1 = varargin{1};
        eps2 = varargin{2};
        eps3 = varargin{3};
        eps4 = varargin{4};
        nn = varargin{5};
        [~, Dmat] = selecpts(z, eps1, eps2, eps3, eps4); % select pts
        Aest = BSS_Moreau(Dmat,N,nn); % BSS dessus
        Best = inv(Aest);
        vectau_init = 1;
        pileB_init = Best;
     case 'truemat'
        T = size(z,2);
        pileB = varargin{1};
        vectau_init = 1:10:T;
        pileB_init = pileB(:,:,vectau_init);
    otherwise
        error('Initialization method unavailable. Choose ''truemat'', ''sobi'' or ''psobi''.')
end