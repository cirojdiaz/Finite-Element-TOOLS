%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Code to compute the stifness matrix in a mixed Problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function A = FEM_Matrix_QQ(LIM , GFX , GFY , GFXT , GFYT , TX , TY , K , NF , N , M, R, Ng, dx , dy ,Dim)

K_GX    = @(x,y) kron(ones(NF),K{1,1}(x,y)).*kron(ones(1,NF),GFX(TX(x),TY(y))) + ...    % X-componente   (K*Grad_Phi)_x = (K11 * Phi_x  +  K12 * Phi_y)  
              kron(ones(NF),K{1,2}(x,y)).*kron(ones(1,NF),GFY(TX(x),TY(y))); 
K_GY    = @(x,y) kron(ones(NF),K{2,1}(x,y)).*kron(ones(1,NF),GFX(TX(x),TY(y))) + ...    % Y-componente   (K*Grad_Phi)_y = (K21 * Phi_x  +  K22 * Phi_y)  
              kron(ones(NF),K{2,2}(x,y)).*kron(ones(1,NF),GFY(TX(x),TY(y)));               
K_G     = @(x,y) kron(ones(NF,1),GFXT(TX(x),TY(y))) .* K_GX(x,y) + ...                % Construindo a Função (tensor)   Grad_Phi_i * (K*Grad_Phi_j)
              kron(ones(NF,1),GFYT(TX(x),TY(y))) .* K_GY(x,y);
          
APT_A   = IntGauss(K_G, LIM(1),LIM(2),LIM(3),LIM(4),N,M,'S',Ng) * (4/(dx*dy));                 % Matriz das aportações de A .................. Integrais dos produtos (Grad_Phi_i) * (K*Grad_Phi_j)

%%%%%%%%%%%%%%%%%%%%%% Construindo os Indices %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I   = transpose(reshape((1:N*M)-1,N,M));                                        % Matris dos N^2 indices ordenados en forma crescente
IA  = kron(ones(NF),I * NF^2);                                                  % Bloques (4x4) dos incrementos de cada indice
ISA = kron(reshape(1:NF^2,NF,NF),ones(M,N));                                    % Incrementos em blocos 4x4 
RX  = kron(  R  ,ones(1,NF));                                                   % X-coneção repetida para cada retangulo
RY  = kron(R(:)',ones(NF,1));                                                   % Y-coneção repetida para cada retangulo
%%%%%%%%%%%%%%%%%%%%%% Construindo a Matriz   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A       = sparse(RX(IA+ISA) , RY(IA+ISA) , APT_A , Dim , Dim); 

end
