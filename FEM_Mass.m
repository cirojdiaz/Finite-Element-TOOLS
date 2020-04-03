%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Computing Mass Matrix in 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = FEM_Mass(LIM ,  TOOLS)


M_Tensor  = @(x,y) kron(ones(1,TOOLS.NF),TOOLS.F(TOOLS.TX(x),TOOLS.TY(y))) ...
               .*  kron(ones(TOOLS.NF,1),TOOLS.FT(TOOLS.TX(x),TOOLS.TY(y)));        % M_F ---> Função (tensor)   Phi_i * Phi_j)

APT_Loc   = FEM_IntGauss(M_Tensor, LIM(1) , LIM(1)+TOOLS.dx , LIM(3) , LIM(3) +  TOOLS.dy , 1 , 1 ,'S',TOOLS.Ng); 
APT_M     = kron(ones(1,TOOLS.M * TOOLS.N),APT_Loc);
%%%%%%%%%%%%%%%%%%%%%% Construindo os Indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IX  =  kron(TOOLS.R , ones(1,TOOLS.NF));
IY  =  kron(TOOLS.R(:)' , ones(TOOLS.NF,1));
%%%%%%%%%%%%%%%%%%%%%% Construindo a Matriz   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M       = sparse(IX,IY, APT_M, TOOLS.Dim , TOOLS.Dim); 

% spy(M);
% full(M)

end