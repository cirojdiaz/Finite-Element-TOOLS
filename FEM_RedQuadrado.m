%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          Funccao para calcular a CONECTIVIDADE de FEM e FV    
%%%
%%%  R --> Matriz de conetividade  ||  pos  --> posiccao dos bordos malha-primal
%%%  Z --> (X,Y) pontos da malha   ||  posV --> posiccao dos bordos malha-dual
%%%  Nx --> # de pontos por arista ||  Nv --> # de volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos,Z,R,Nx,Nr,NF,X_nodes,Y_nodes,pos_Herm,X_Herm,Y_Herm,Inner_Herm] = FEM_RedQuadrado(LIM,n,m,Q_Type)

if nargin ==1
    Q_Type=1;                                                               % Ordem de Q_FEM   (Q1 , Q2 , Q3 ....) 
end   

if Q_Type == 1 
    ord = [1;2;3;4];
elseif Q_Type == 2
%     ord=[1;7;9;3;4;8;6;2;5];
    ord = (1:9)';
elseif Q_Type == 3   
    ord = (1:16)';
elseif Q_Type == 4   
    ord = [ 1;6;11;16;21;2;7;12;17;22;3;8;13;18;23;4;9;14;19;24;5;10;15;20;25]; 
%     ord = 1:25;
end

NF   = (Q_Type+1)^2;
Nx   = 2^n * Q_Type + 1;                                                    % Nx --> # de pontos por arista   
Ny   = 2^m * Q_Type + 1;                                                    % Nx --> # de pontos por arista   
Npx   = 2^n; 
Npy   = 2^m; 
Nr   = Npx * Npy;
X    = linspace(LIM(1),LIM(2),Nx); 
Y    = linspace(LIM(3),LIM(4),Ny);
[X,Y]= meshgrid(X,Y);                                                       % Malha Nx X Nx em [0,1]x[0,1]
Z    = [reshape(Y,1,Nx*Ny);reshape(X,1,Nx*Ny)];                               % Calculando las cordenadas de los nodos
Num  = (1:Nx*Ny); 
Num = reshape(Num,[Nx,Ny])'; 
pos  = {Num(1:Ny,1);Num(Ny,1:Nx)';Num(Ny:-1:1,Nx);Num(1,Nx:-1:1)'};         % Almazenando posiciones de la matris corrspondientes a los Bordes
X_nodes = {X(:,1) ; X(end,:)' ; X(end:-1:1,end) ; X(1,end:-1:1)'}; 
Y_nodes = {Y(:,1) ; Y(end,:)' ; Y(end:-1:1,end) ; Y(1,end:-1:1)'}; 
%%%% Special case for Hermite Polynomials
if Q_Type == 3
    pos_Herm = {Num(1:3:Ny,1);Num(Ny,1:3:Nx)';Num(Ny:-3:1,Nx);Num(1,Nx:-3:1)'};
    Inner_Herm = Num';  
    Inner_Herm(unique(cell2mat(pos_Herm))) = [];
    X_Herm = {X(1:3:end,1) ; X(end,1:3:end)' ; X(end:-3:1,end) ; X(1,end:-3:1)'}; 
    Y_Herm = {Y(1:3:end,1) ; Y(end,1:3:end)' ; Y(end:-3:1,end) ; Y(1,end:-3:1)'}; 
end
%%%%%%%%%%%%%%%%%%  Imprimir a malha %%%%%%%%%%%%%%%%%%%%%%%%
% text(X(:),Y(:),num2str(Num(:)),'color','b','fontsize',10);
% axis([-0.1,1.1,-0.1,1.1]); 
% set(gca,'xtick',[],'ytick',[]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R  = zeros(NF,Nr);                                                        
RP = kron(ones(Q_Type+1,1),(0:Q_Type)) +...                                 % Ordem na que soman os elementos de cada quadrado
     Nx*kron((0:Q_Type)',ones(1,Q_Type+1));                                 % Q1: RP = kron([1;1],[0,1])  +  Nx*kron([0;1],[1,1]);
                                                                            % Q2: RP = kron([1;1;1],[0,1,2]) + Nx*kron([0;1;2],[1,1,1]);  
for j = 1:Npy
    for i = 1:Npx
        N       = Nx*(j-1)+i;
        Rtemp   = RP' + (Q_Type)*N - (Q_Type-1);                                              
        Rcol    = Rtemp(ord);                                               % Organizando as colunas do retangolo pelo ordem das funções bases locais 
        R(:,Npx*(j-1)+i)  = Rcol;
    end
end

% disp(num2str(R));

end















