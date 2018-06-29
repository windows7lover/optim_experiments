function out = solve_lasso_cvx(A,b,lambda)

[n,d] = size(A);

cvx_begin
    cvx_solver sdpt3
    variables x(d,1) t(d,1)
    minimize 0.5*sum_square(A*x-b) + lambda * norm(x,1)

cvx_end
out = x;