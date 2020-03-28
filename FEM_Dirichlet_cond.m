%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Code to imposse Dirichlet Boundary conditions to an Element S
%               (It could either be a matrix or a vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S,Rows] = FEM_Dirichlet_cond(S,pos,TOOLS, G , Rows)

N = size(S,1);
M = size(S,2);

%%% Simple Homogeneous %%%
if nargin == 2
    pos = unique(cell2mat(pos));
    if N==1 || M == 1
        S(pos) = [];
    else
        S(pos,:) = [];
        S(:,pos) = [];
    end
else
    %%%%  General Dirichlet
    [pos,ind] = unique(cell2mat(pos));
    if N==1 || M == 1
        x_nodes = cell2mat(TOOLS.X_nodes);
        y_nodes = cell2mat(TOOLS.Y_nodes);
        S = S - Rows * G(x_nodes(ind) , y_nodes(ind));
        S(pos) = G(x_nodes(ind) , y_nodes(ind));
%         S(pos) = [];
    else
        Rows = S(:,pos);
%         S(pos,:) = [];
%         S(:,pos) = [];
        S(pos,:) = 0;
        S(:,pos) = 0;
        S(sub2ind([N,M],pos,pos)) = 1;
    end
end