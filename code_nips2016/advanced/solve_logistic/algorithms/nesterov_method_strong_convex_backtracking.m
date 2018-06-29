function [xout,param] = nesterov_method_strong_convex_backtracking(finfo,xHist,k,param)

% L = finfo.L;
mu = finfo.mu;
f = finfo.f;
fp = finfo.fp;
yold = xHist(:,k-1);


if(k == 2)
    param.xold = yold;
    param.L = finfo.L;
end

L = param.L;
Q = (sqrt(L)-sqrt(mu))/(sqrt(L)+sqrt(mu));

xold = param.xold;
fpyold = fp(yold);
fyold = f(yold);
while(true)
    xplus = yold - fpyold/L;
    if( f(xplus) > fyold + fpyold'*(xplus-yold) + (L/2)*norm(xplus-yold)^2 ) % Problem: constant too small
        L = L*2;
    else
        break;
    end
end

yplus = xplus + Q*(xplus-xold);
param.xold = xplus;

xout = yplus;
param.L = L/2;