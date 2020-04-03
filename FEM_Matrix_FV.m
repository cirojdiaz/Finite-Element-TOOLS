%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Computing the Finite Volumes Matrix 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AR = FEM_Matrix_FV(LIM , GFX , GFY , TX , TY ,  Nvx , Nvy , NF , N , M, R, K , Ng, Dim)

VS    = ones(NF,1); 

Gx    = @(x,y) kron(VS,K{1,1}(x,y)).*kron(ones(1,T),GFX(TX(x),TY(y))) + ...      % X-componente   (K*Grad_Phi)_x = (K11 * Phi_x  +  K12 * Phi_y)  
               kron(VS,K{1,2}(x,y)).*kron(ones(1,T),GFY(TX(x),TY(y))); 
Gy    = @(x,y) kron(VS,K{2,1}(x,y)).*kron(ones(1,T),GFX(TX(x),TY(y))) + ...      % Y-componente   (K*Grad_Phi)_y = (K21 * Phi_x  +  K22 * Phi_y)  
                kron(VS,K{2,2}(x,y)).*kron(ones(1,T),GFY(TX(x),TY(y))); 


GL        = @(x,y) [Gx(x,y) , Gy(x,y)];
K_GX_R    = @(x,y) -(2/dx) * removerows(transpose(GL(x,y)),size(x,2)+1:2*size(x,2))';
K_GY_R    = @(x,y) -(2/dy) * removerows(transpose(GL(x,y)),1          : size(x,2) )';

APT_ARX = IntGauss(K_GX_R,LIM(1),LIM(2),LIM(3),LIM(4),2*N,2*M,'L',Ng);                       % Matriz das aportações de AR .................. Integrais de linha de (K*Grad_Phi_j) nas arestas dos volumes
APT_ARY = IntGauss(K_GY_R,LIM(1),LIM(2),LIM(3),LIM(4),2*N,2*M,'L',Ng);                       % Matriz das aportações de AR .................. Integrais de linha de (K*Grad_Phi_j) nas arestas dos volumes 

%%%%%%%%%%%%%%%%%%%%%% Construindo Indices AR e g  %%%%%%%%%%%%%%%%%%%%%%%%
Tx       =  [2:2:2*N,(2:2:2*N) + 2*N-1];                                     % Indices pares, usamos para reter so as colunas pares de ARX e ARY 
Ty       =  kron(ones(NF,1),[2:2:2*M,(2:2:2*M) + 2*M-1]') +...
              (4*M)*kron((0:NF-1)',ones(2*M,1));                             % Indices pares, usamos para reter so as colunas pares de ARX e ARY                  
ARX     =   -APT_ARX{2}(:,Tx);                                                % Tirando as colunas imapars em APT_ARX em APT_ARY 
ARY     =   APT_ARY{1}(Ty,:);                                                % Tirando as colunas imapars em APT_ARY 
ARX(:,end/2+1:end) = -ARX(:,end/2+1:end);                                    % Trocando sinal de acordo ao sentido da integral en Dx
Ind_Sgn = kron(ones(1,NF),1:M) + 2*M*kron((0:NF-1),ones(1,M));
ARY(Ind_Sgn,:) = -ARY(Ind_Sgn,:);                                            % Trocando sinal de acordo ao sentido da integral en Dy

VX = kron(ones(NF,1),kron([(1:N)+1,1:N],ones(2*M,1)) +...
                     kron(kron(ones(M,1),[0;N+1]) +...
                      (N+1)*kron((0:M-1)',[1;1]) , ones(1,2*N)));
VY = kron(ones(NF,1),kron(ones(1,2*N), kron(ones(2,1),(1:N+1:N*M)') +...
                    (N+1)*kron([1;0],ones(M,1))) + ...
                    kron(ones(2*M,1),intercale(0:N-1,1:N)));
                
BLX  = kron(ones(M,1),kron(ones(2),1:N)) +  N*kron((0:M-1)',ones(2,2*N));        
t1x   = kron((1:NF)',ones(size(BLX)));
t2x   = kron(ones(NF,1),BLX);
RRX = R(sub2ind([NF,M*N],t1x,t2x));

BLY  = kron(ones(2,1),kron((reshape((1:N*M),N,M))',ones(1,2)));        
t1y   = kron((1:NF)',ones(size(BLY)));
t2y   = kron(ones(NF,1),BLY);
RRY = R(sub2ind([NF,M*N],t1y,t2y));

%%%%%%%%%%%%%%%%%%%%%% Construindo a Matriz   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AR      = sparse([VX,VY], [RRX,RRY], [ARX,ARY] , Nvx*Nvy , Dim );                   % Construindo a matriz de restrições. 

end

function G = K_Grad(K,x,y,Type,n,m,medio)
global NF GFX GFY LIM
x1 = LIM(1); x2 = LIM(2); y1 = LIM(3); y2 = LIM(4);
N   = 2^n;                        M   = 2^m;                                 % Numero de retángulos ....  N*M
dx  = (x2-x1)/N;                  dy  = (y2-y1)/M;                           % Norma da partição em X e Y;
TX  = @(x)  2*((x-x1)/dx - floor((x-x1)/dx))-1;                              % X-Transformação que leva no retangulo padrão [-1,1]x[-1,1]
TY  = @(y)  2*((y-y1)/dy - floor((y-y1)/dy))-1;                              % Y-Transformação que leva no retangulo padrão [-1,1]x[-1,1]                 
gfx   = GFX(TX(x),TY(y));       gfy   = GFY(TX(x),TY(y));
T     = (Type=='S')*NF + (Type == 'L')*1;
VS    = ones(NF,T); 
if strcmp(medio,'Anisotropico')
Gx    = kron(VS,K{1,1}(x,y)).*kron(ones(1,T),gfx) + ...      % X-componente   (K*Grad_Phi)_x = (K11 * Phi_x  +  K12 * Phi_y)  
        kron(VS,K{1,2}(x,y)).*kron(ones(1,T),gfy); 
Gy    = kron(VS,K{2,1}(x,y)).*kron(ones(1,T),gfx) + ...      % Y-componente   (K*Grad_Phi)_y = (K21 * Phi_x  +  K22 * Phi_y)  
        kron(VS,K{2,2}(x,y)).*kron(ones(1,T),gfy); 
G     = [Gx,Gy];
elseif strcmp(medio,'Heterogenous')
    Gx    = kron(VS,K).*kron(ones(1,T),gfx); 
    Gy    = kron(VS,K).*kron(ones(1,T),gfy);
    G     = [Gx,Gy];
else
    G = NaN;
end
end