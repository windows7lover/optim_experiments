function X = linear_fixed_point(A,xstar,x0,k)
% Does k iteration x=Ax+b and return the k+1 points [x0, x1,..., xk].
% b is computed s.t. xstar = Axstar-b

x = x0;

X = [x, zeros(length(x),k)];

for i=1:k
    x = A*(x-xstar)+xstar;
    X(:,i+1) = x;
end