%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Computing the Drift Matrix for parabolic equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function D = FEM_Drift(LIM , TOOLS , prmt)

if nargin==2
    prmt = struct;
    prmt.b = [1,1];
end

D_Tensor     = @(x,y) prmt.b(1)*(2/TOOLS.dx)*kron(ones(TOOLS.NF,1),TOOLS.GFXT(TOOLS.TX(x),TOOLS.TY(y)))...
                                          .* kron(ones(1,TOOLS.NF),TOOLS.F(TOOLS.TX(x),TOOLS.TY(y)))  + ...                % Construindo a Função (tensor)   Grad_Phi_i * (K*Grad_Phi_j)
                      prmt.b(2)*(2/TOOLS.dy)*kron(ones(TOOLS.NF,1),TOOLS.GFYT(TOOLS.TX(x),TOOLS.TY(y))) ...
                                          .* kron(ones(1,TOOLS.NF),TOOLS.F(TOOLS.TX(x),TOOLS.TY(y)));

% D_Tensor     = @(x,y) prmt.b(1)*(2/TOOLS.dx)*kron(ones(1,TOOLS.NF),TOOLS.F(TOOLS.TX(x),TOOLS.TY(y)))...
%                                           .* kron(ones(TOOLS.NF,1),TOOLS.GFXT(TOOLS.TX(x),TOOLS.TY(y)))  + ...                % Construindo a Função (tensor)   Grad_Phi_i * (K*Grad_Phi_j)
%                       prmt.b(2)*(2/TOOLS.dy)*kron(ones(1,TOOLS.NF),TOOLS.F(TOOLS.TX(x),TOOLS.TY(y))) ...
%                                           .* kron(ones(TOOLS.NF,1),TOOLS.GFYT(TOOLS.TX(x),TOOLS.TY(y)));
                  
APT_Loc   = FEM_IntGauss(D_Tensor, LIM(1) , LIM(1)+TOOLS.dx , LIM(3) , LIM(3)+TOOLS.dy ,1,1,'S',TOOLS.Ng); 
APT_D     = kron(ones(1,TOOLS.M * TOOLS.N),APT_Loc);
%%%%%%%%%%%%%%%%%%%%%% Construindo os Indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IX  =  kron(TOOLS.R,ones(1,TOOLS.NF));
IY  =  kron(TOOLS.R(:)',ones(TOOLS.NF,1));
%%%%%%%%%%%%%%%%%%%%%% Construindo a Matriz   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D       = sparse(IX,IY, APT_D, TOOLS.Dim , TOOLS.Dim); 

% spy(D);
% full(D)

end