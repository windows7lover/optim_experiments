function [A, b, xsol] = gen_problem_Lasso_random(n,d,s_approx,cond_number,norm_xstar,lambda)

if(s_approx >= d)
    error('s_approx should be < d')
end

A=randn(n,d);
[U,S,V]=svd(A);
S(S~=0)=1+rand(1,min(d,n))*cond_number;
A=U*S*V';

x = rand(d,1)-0.5;
x = norm_xstar*x/norm(x);
b = A*x;

xsol = solve_lasso_cvx(A,b,lambda);