function [ f_seq ] = nesterov( x0, finfo, kmax )

f = finfo.f;
d = finfo.n;
fstar = finfo.fstar;

x = x0;
y = x0;

kappa = finfo.mu/finfo.L;
beta = (1-sqrt(kappa))/(1+sqrt(kappa));

f_seq = zeros(1,kmax+1);
f_seq(1,1) = f(x0) - fstar;

Lk = finfo.L;

for i = 2:kmax+1
    % xnew = y-finfo.fp(y)/finfo.L;
    xnew = damped_gradient(y, finfo, Lk);
    y = (1+beta) * xnew - beta *x;
    x = xnew;
    f_seq(1,i) = f(x)-fstar;
end

end

