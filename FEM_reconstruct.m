%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Code to reconstruct the solution in a mesh [X,Y] from
%   the coeficients of FEM Qr ... r =1,2,3 and Hermite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SOL = FEM_reconstruct(Omega,Omega_e,C,X,Y,TOOLS)

e = 10^-10;
X = X - (rem(X,TOOLS.dx) == 0).*(X~=Omega_e(1))*e;      % Keeping function TX good evaluated at the end of each element
Y = Y - (rem(Y,TOOLS.dx) == 0).*(Y~=Omega_e(3))*e;      % Keeping function TY good evaluated at the end of each element

% N_nodes = sqrt(numel(C));
% C = padarray(reshape(C , N_nodes , N_nodes),[1,1],0,'both');

[n,m] = size(X);

XLoc = floor((X-Omega(1))./TOOLS.dx) +1;
IXR = XLoc.*(XLoc <= TOOLS.N) + TOOLS.N.*(XLoc > TOOLS.N);
YLoc = floor((Y-Omega(3))./TOOLS.dy) +1;
IYR = YLoc.*(YLoc <= TOOLS.M) + TOOLS.M.*(YLoc > TOOLS.M);

R_Elmnt  = kron(ones(TOOLS.NF,1) , IXR + (IYR-1) * TOOLS.M);
R_Fun  = kron((1:TOOLS.NF)',ones(size(X)));

% tmp1 = sub2ind([TOOLS.NF , TOOLS.Dim],R_Fun,R_Elmnt);
% tmp3 = TOOLS.R(tmp1);

IND = TOOLS.R(sub2ind([TOOLS.NF , TOOLS.N * TOOLS.M],R_Fun,R_Elmnt));
SUM = kron(ones(1,TOOLS.NF) , speye(n,m));

% tmp2 = C(IND);
% tmp3 = TOOLS.F( TOOLS.TX(X) , TOOLS.TY(Y));
% figure; surf(tmp2);

SOL = SUM * ( C(IND) .* TOOLS.F( TOOLS.TX(X) , TOOLS.TY(Y)) );

% figure; surf(TOOLS.F( TOOLS.TX(X) , TOOLS.TY(Y)));

end













    