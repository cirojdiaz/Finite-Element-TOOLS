%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Codigo para el Calculo de las derivadas simbolicas de una funcion
%                           FUN(T,X,Y)   
%            REMPLAZANDO parametros por su valores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [QT,QX,QY,QXX,QYY,Q_Lap] = FEM_Derivatives(pname,pval,Fun,VAR)

n_pmt = numel(pname);
values = cell(min(1,n_pmt),n_pmt);

for j = 1:n_pmt
    values{j} = num2str(pval{j});
end

if nargin==3
   VAR = '@(x,y) '; 
end

V = str2sym(func2str(Fun));
syms t x y

QT      = str2func([VAR,replace(extractAfter(func2str(matlabFunction(simplify(diff(V,t)))),')'),pname,values)]);
QX      = str2func([VAR,replace(extractAfter(func2str(matlabFunction(simplify(diff(V,x)))),')'),pname,values)]);
QY      = str2func([VAR,replace(extractAfter(func2str(matlabFunction(simplify(diff(V,y)))),')'),pname,values)]);
QXX     = str2func([VAR,replace(extractAfter(func2str(matlabFunction(simplify(diff(diff(V,x),x)))),')'),pname,values)]);
QYY     = str2func([VAR,replace(extractAfter(func2str(matlabFunction(simplify(diff(diff(V,y),y)))),')'),pname,values)]);

Q_Lap       = @(x,y) -QXX(x,y) - QYY(x,y);

end