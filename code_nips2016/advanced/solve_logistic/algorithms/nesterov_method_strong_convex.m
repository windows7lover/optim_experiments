function [xout,param] = nesterov_method_strong_convex(finfo,xHist,k,param)

L = finfo.L;
mu = finfo.mu;
fp = finfo.fp;
yold = xHist(:,k-1);

beta = (sqrt(L)-sqrt(mu))/(sqrt(L)+sqrt(mu));

if(k == 2)
    param.xold = yold;
end

xold = param.xold;
xplus = yold - fp(yold)/L;
yplus = xplus + beta*(xplus-xold);
param.xold = xplus;

xout = yplus;