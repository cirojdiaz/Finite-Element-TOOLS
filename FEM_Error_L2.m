%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Function to compute L2 errors between Coeficient Solutions
%    C, Sol in R2 ... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = FEM_Error_L2(C , Sol , LIM , LIM_E, TOOLS, TOOLS_last)

%%%%%%  Meshing the coeficient  %%%%%%
      
%%%%% Find Nodes in the efective area %%%%%
x      = linspace(LIM(1),LIM(2),TOOLS.Nx);
[~,ix] = find( (x >= LIM_E(1)) .* (x<=LIM_E(2)) );
y      = linspace(LIM(3),LIM(4),TOOLS.Ny);
[~,iy] = find( (y >= LIM_E(1)) .* ( y<=LIM_E(2) ));

[X,Y] = meshgrid( x , y );
if isa(Sol,'function_handle')
    S     = Sol(X,Y);
elseif isa(Sol,'double')
    S  = FEM_reconstruct(LIM,LIM,Sol',X,Y,TOOLS_last);
end
 
% figure; surf(S);
% figure; surf(C);
C(ix,:) = 0; C(:,iy) = 0; 
S(ix,:) = 0; S(:,iy) = 0;

E = sqrt( ((C(:) - S(:))' * TOOLS.Mass * (C(:) - S(:)) )  / ( S(:)' * TOOLS.Mass * S(:) ) );




end