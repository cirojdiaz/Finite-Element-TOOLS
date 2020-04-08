function b = FEM_Binomial_Coef(n,b)

if n == 0
    b = 1;
    return
else 
    if nargin == 1
        b = 1;
    end
    b = [FEM_Binomial_Coef(n-1,b),0] + [0,FEM_Binomial_Coef(n-1,b)];
end

end