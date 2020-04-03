%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Computing the line integrals from part integration 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = FEM_Line_Integral(LIM , TOOLS , prmt)

if nargin==2
    prmt = struct;
    prmt.sigma = sqrt(2)*[1,0;0,1];
end

e = 10^-15;

L_TensorX  = @(x,y) prmt.sigma(2,2)^2/2 * (2/TOOLS.dx) * kron(ones(TOOLS.NF,1),TOOLS.FT(TOOLS.TX(x),TOOLS.TY(y))) ...
                                                      .* kron(ones(1,TOOLS.NF),TOOLS.GFY(TOOLS.TX(x),TOOLS.TY(y)));  
            
L_TensorY  = @(x,y) prmt.sigma(1,1)^2/2 * (2/TOOLS.dy) * kron(ones(TOOLS.NF,1),TOOLS.FT(TOOLS.TX(x),TOOLS.TY(y))) ...
                                                      .* kron(ones(1,TOOLS.NF),TOOLS.GFX(TOOLS.TX(x),TOOLS.TY(y)));

[~,~,~,~,px,py,wx,wy] = FEM_GaussMat(LIM(1) , LIM(1)+TOOLS.dx , LIM(3), LIM(3) + TOOLS.dy , TOOLS.Ng);

P1 = - L_TensorY(px                  ,  LIM(3)            );
P2 =   L_TensorX(LIM(1)+ TOOLS.dx - e,    py              );        P2 = P2 .* (abs(P2)>10^-10);
P3 =   L_TensorY(py                  ,LIM(3)+TOOLS.dy - e );        P3 = P3 .* (abs(P3)>10^-10);
P4 = - L_TensorX(LIM(1)              ,    py              );

APT_Loc1  =  P1 * kron(speye(TOOLS.NF),wx');
APT_Loc2  =  P2 * kron(speye(TOOLS.NF),wy');
APT_Loc3  =  P3 * kron(speye(TOOLS.NF),wx');
APT_Loc4  =  P4 * kron(speye(TOOLS.NF),wy');

% Intersections = ((APT_Loc1~=0) + (APT_Loc2~=0) + (APT_Loc3~=0) + (APT_Loc4~=0))==2;

APT_Loc = (APT_Loc1 + APT_Loc2 + APT_Loc3 + APT_Loc4);%.*(Intersections*(-1/2) + 1);

% APT_M     = kron(ones(1,2 * (TOOLS.M + TOOLS.N)),APT_Loc);
APT_L     = kron(ones(1,TOOLS.M * TOOLS.N),APT_Loc);

%%%%%%%%%%%%%%%%%%%%%% Construindo os Indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOOLS.R(sqrt(TOOLS.NF)+1:end,:) = [];
IX  =  kron(TOOLS.R , ones(1,TOOLS.NF));
IY  =  kron(TOOLS.R(:)' , ones(TOOLS.NF,1));
%%%%%%%%%%%%%%%%%%%%%% Construindo a Matriz   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pos      = unique(cell2mat(TOOLS.pos));
[posX,posY] = meshgrid(pos);
Border   = sparse(posX , posY , 1 , TOOLS.Dim , TOOLS.Dim);
L       = sparse(IX,IY, APT_L, TOOLS.Dim , TOOLS.Dim) .* Border; 
% spy(L)
end