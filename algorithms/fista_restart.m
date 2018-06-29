function [xout,param] = fista_restart(finfo,x,k,param)

L = finfo.L;
f = finfo.f;
fp = finfo.fp;

theta = @(ite) 2/(ite+3);

if(k == 2)
    param.xold = x;
    theta = @(ite) 1;
    param.fold = inf;
end

xold = param.xold;

h = 1/L;
coef = theta(k)*( (1/theta(k-1)) -1 );
y = x+coef*(x-xold);
yplus = y - h*fp(y);
xout = finfo.proxoperator.prox(yplus,h);

param.xold = x;

% restart test
fout = f(xout) + finfo.proxoperator.f(xout);
if(fout>param.fold)
    param.xold = xout;
end
param.fold = fout;