function [heapB_init, vectau_init] = JEFASBSSinit(z, init_meth, varargin)
% JEFASBSSinit Choose BSS method for initialization of JEFAS-BSS
% usage:	[heapB_init, vectau_init] = JEFASBSSinit(z, init_meth, varargin)
%
% Input:
%   z: observed signals
%   inith_meth: initialization method. MUST be 'sobi', 'psobi', 'QTFBSS', or 'ground_truth'
%   varargin: 
%       if inith_meth = 'sobi': varargin unused
%       if inith_meth = 'psobi': varargin=Dt subsampling time for unmixing matrices estimations
%       if inith_meth = 'QTFBSS': varargin={eps3,eps4,nn} see BSS_QTF.m 
%       if inith_meth = 'ground_truth': varargin=heapB Ground truth unmixing matrices
%
% Output:
%   heapB_init: initial unmixing matrices (third dimension for time)
%   vectau_init: time where the initial unmixing matrices are given

% Copyright (C) 2018 Adrien MEYNARD
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% 
% Author: Adrien MEYNARD
% Created: 2018-10-30

switch init_meth
    case 'sobi'
        N = size(z,1);
        [Asobi,~] = sobi(z,N,1000);
        B0 = inv(Asobi);
        B0 = B0./sqrt(sum(abs(B0).^2,2));
        vectau_init = 1;
        heapB_init = B0;
    case 'psobi'
        T = size(z,2);
        Dt = varargin{1};
        vectauns = 1:Dt:T;
        [~, heap_Bsobins] = psobi(z,vectauns);
        vectau_init = vectauns;
        heapB_init = heap_Bsobins;
    case 'QTFBSS'
        N = size(z,1);
        eps3 = varargin{1};
        eps4 = varargin{2};
        nn = varargin{2};
        [~, Dmat] = selecpts(z, 1, 0, eps3, eps4); % select pts
        Aest = BSS_Moreau(Dmat,N,nn); % BSS dessus
        Best = inv(Aest);
        vectau_init = 1;
        heapB_init = Best;
     case 'ground_truth'
        T = size(z,2);
        heapB = varargin{1};
        vectau_init = 1:10:T;
        heapB_init = heapB(:,:,vectau_init);
    otherwise
        error('Initialization method unavailable. Choose ''ground_truth'', ''sobi'',''psobi'', or ''QTFBSS''.')
end