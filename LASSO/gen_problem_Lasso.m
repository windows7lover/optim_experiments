function [A, b, xsol] = gen_problem_Lasso(rho, lambda, s, d, n)
% This generates A, b, x_sol, sol such that x_sol is solution of
% min_x 0.5*||Ax-b||_2^2 + lambda ||x||_1
% and sol is the optimal value
%
% This is taken form gradient methods for minimizing composite function of
% nesterov page 25-26
%
% n,d dimensions of the problem
% s sparisty of the built solution
% rho controls the magnitude of x_sol (take rho = 1 for example)
%
% n >= floor(s*log(d)) to get a solution

B = 2*rand(n,d)-1;

v = rand(n,1);

y = v/norm(v,2);

[u,idx] = sort(B'*y,'descend');

A = zeros(n,d);
for i  = 1:d
    if i<=s
        A(:,i) = B(:,idx(i))/abs(u(i));
    else
        if i>s
            if abs(u(i))<=0.1
                A(:,i)= B(:,idx(i));
            else
                A(:,i) = rand*B(:,idx(i))/abs(u(i));
            end
        end
    end
end

A = A*lambda;

v  = A'*y;
xsol = zeros(d,1);
for i = 1:s
    xsol(i) = (rho/sqrt(s))*rand*sign(v(i));
end

b = y+A*xsol;
end