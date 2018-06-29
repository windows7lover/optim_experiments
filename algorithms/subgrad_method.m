function [xplus, param] = subgrad_method(fun,xHist,k,param)

gx = fun.fp;
x = xHist(:,k-1);

if(isfield(param,'nite'))
    param.nite = param.nite+1;
    h = 1/param.nite;
else
    h = 1/sqrt(k+1);
end

xplus = x-h*gx(x);