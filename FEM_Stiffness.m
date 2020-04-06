%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Code to compute the Stifness matrix of FEM using bassis
%               Qr ... r =1,2,3 and Hermite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = FEM_Stiffness(LIM , TOOLS , prmt)

if nargin==2
    prmt = struct;
    prmt.sigma = sqrt(2) * [1,0;0,1];
end
%%%%%%%%%% Construindo a Função Tensor %%%%%%%%%%%%%
A_Tensor     = @(x,y) prmt.sigma(1,1)^2/2 * (4/TOOLS.dx^2)*kron(ones(TOOLS.NF,1),TOOLS.GFXT(TOOLS.TX(x),TOOLS.TY(y))) ...
                                                  .* kron(ones(1,TOOLS.NF),TOOLS.GFX(TOOLS.TX(x),TOOLS.TY(y))) + ...      
                      prmt.sigma(2,2)^2/2 * (4/TOOLS.dy^2)*kron(ones(TOOLS.NF,1),TOOLS.GFYT(TOOLS.TX(x),TOOLS.TY(y))) ...
                                                  .* kron(ones(1,TOOLS.NF),TOOLS.GFY(TOOLS.TX(x),TOOLS.TY(y)));
% A_Tensor     = @(x,y) prmt.sigma(1,1)^2/2 * (4/TOOLS.dx^2)*kron(ones(1,TOOLS.NF),TOOLS.GFX(TOOLS.TX(x),TOOLS.TY(y))) ...
%                                                   .* kron(ones(TOOLS.NF,1),TOOLS.GFXT(TOOLS.TX(x),TOOLS.TY(y))) + ...      
%                       prmt.sigma(2,2)^2/2 * (4/TOOLS.dy^2)*kron(ones(1,TOOLS.NF),TOOLS.GFY(TOOLS.TX(x),TOOLS.TY(y))) ...
%                                                   .* kron(ones(TOOLS.NF,1),TOOLS.GFYT(TOOLS.TX(x),TOOLS.TY(y)));

APT_Loc   = FEM_IntGauss(A_Tensor, LIM(1) , LIM(1)+TOOLS.dx , LIM(3) , LIM(3)+TOOLS.dy , 1 , 1 , 'S' , TOOLS.Ng) ; 
APT_A     = kron(ones(1,TOOLS.M * TOOLS.N),APT_Loc);
%%%%%%%%%%%%%%%%%%%%%% Construindo os Indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IX  =  kron(TOOLS.R,ones(1,TOOLS.NF));
IY  =  kron(TOOLS.R(:)',ones(TOOLS.NF,1));
%%%%%%%%%%%%%%%%%%%%%% Construindo a Matriz   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A       = sparse(IX,IY, APT_A, TOOLS.Dim , TOOLS.Dim); 

% spy(A);
% full(A)

end