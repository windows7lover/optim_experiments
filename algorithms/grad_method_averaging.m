function [xplus, param] = grad_method_averaging(fun,x,k,param)

L = fun.L;
gx = fun.fp;

if(k == 2)
    param.xgrad= x;
end

xgrad_old = param.xgrad;

xgrad = xgrad_old-(1/L)*gx(xgrad_old);
param.xgrad = xgrad;

xplus = (xgrad + x * k)/(k+1);
