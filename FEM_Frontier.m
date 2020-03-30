%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Vector containing the Frontier values of function G
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = FEM_Frontier(G,TOOLS)

[pos,ind] = unique(cell2mat(TOOLS.pos));
x_nodes = cell2mat(TOOLS.X_nodes);
y_nodes = cell2mat(TOOLS.Y_nodes);

b = sparse(pos , ones(numel(ind),1), G(x_nodes(ind) , y_nodes(ind)) , TOOLS.Dim , 1);

end