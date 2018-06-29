function [algo_x,algo_tolf] = cumul_rmpe(X,lambda,finfo)

[d,n] = size(X);

algo_x = zeros(d,n-2);
algo_tolf = zeros(n-2,1);

for i=3:n
    U = diff(X(:,1:i),1,2);
    UU = U'*U;
    k = i-1;
    UU = UU/norm(UU);
    
    z = (UU+lambda*eye(k))\ones(k,1);
    c = z/sum(z);
    algo_x(:,i-2) = X(:,1:(i-1))*c;
    algo_tolf(i-2) = finfo.f(algo_x(:,i-2));
end