%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%  PROGRAMA PARA CALCULAR LAS INTEGRALES  DE %%
%%  UMA FUNÇÃO f NO RETANGULO [x1,x2]X[y1,y2] %%
%%  DIVIDIDO EM N,M BLOCOS                    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function II = FEM_IntGauss(Fun,x1,x2,y1,y2,N,M,Type,Ng)
 
e = 10^(-15);                                                               % Shift paar integrais em interiores e exteriores de retangulos
if nargin == 7;  Type = 'S';  end
if nargin == 5; N = 1; M = 1; end

[DX,DY] = size(Fun(x1,y1));                                                   % Calculando a dimenção do espaço de chegada de @f
tx      = (x2-x1)/N;             ty = (y2-y1)/M;                            % Calculando contração do NOVO retangulo padrão 

if strcmp(Type,'S')                                                         %%%%%%%%%  Caso integral de Superficie 
    [px,py,wx,wy] = FEM_GaussMat(x1, x1+tx, y1, y1+ty,Ng);                      % Get gauss points and weights at NOVO Rectangulo padrão
    [n,m]         = meshgrid(1:N, 1:M);                                     % Contructing Gauss index mesh for N*M rectangles
    PX      = kron(ones(size(n)),px) + tx*kron(n-1,ones(size(px)));         % X-Contructing Gauss points mesh along all N*M rectangles
    PY      = kron(ones(size(m)),py) + ty*kron(m-1,ones(size(py)));         % Y-Contructing Gauss points mesh along all N*M rectangles
    WX      = kron(ones(DX,DY),kron(ones(size(n)),wx));                     % X-Contructing Gauss Weights mesh along all N*M rectangles
    WY      = kron(ones(DX,DY),kron(ones(size(m)),wy));                     % Y-Contructing Gauss Weights mesh along all N*M rectangles
    WW = WX.*WY;       clear WX WY
    FF = Fun(PX,PY); 
    PPX     = kron(speye(DX),kron(speye(M),ones(1,Ng)));                    % Matriz para soma por colunas em cada retangulo
    PPY     = kron(speye(DY),kron(speye(N),ones(1,Ng)));                    % Matriz para soma por linhas em cada retangulo  
    II      = PPX*(FF.*WW)*PPY';

elseif strcmp(Type,'L')                                                     %%%%%%%%%  Caso Linea
    [~,~,~,~,px,py,wx,wy] = GaussMat(x1, x1+tx, y1, y1+ty,Ng);                 % Get gauss points and weights at NOVO Rectangulo padrão
    X  = x1:tx:x2;              Y = y1:ty:y2;                               % Contruindo o vetor das Aristas
    X  = [X+e,X-e];  X([N+1,N+2]) = [];                                     % Integrais internas y externas En las X-Aristas
    Y  = [Y+e,Y-e];  Y([N+1,N+2]) = [];                                     % Integrais internas y externas En las Y-Aristas
    PX        = kron(ones(1,N),px) + tx*kron(0:N-1,ones(size(px)));         % N Gauss points on a single X-side
    PY        = kron(ones(1,M),py) + ty*kron(0:M-1,ones(size(py)));         % M Gauss points on a single Y-side
    [PFY,PFX] = meshgrid(PX,Y);                                             % X-Contructing Gauss points mesh along all X-sides
    [PCX,PCY] = meshgrid(PY,X);                                             % X-Contructing Gauss points mesh along all Y-sides
    WX        = kron(ones(DX,DY),kron(ones(2*M,N),wx));                     % Contructing Gauss Weights mesh for points of X-sides
    WY        = kron(ones(DX,DY),kron(ones(2*N,M),wy));                     % Contructing Gauss Weights mesh for points of Y-sides
    PPX       = kron(speye(DY),kron(speye(N),ones(Ng,1)));                   % Matriz para soma por colunas 
    PPY       = kron(speye(DY),kron(speye(M),ones(Ng,1)));                   % Matriz para soma por colunas 
    IIX       = (Fun(PFX,PFY).*WX)*PPX;                                     % Integrando em cada X-side
    IIY       = (Fun(PCX,PCY).*WY)*PPY;                                     % Integrando em cada Y-side
    II        = {IIX,IIY};
end
