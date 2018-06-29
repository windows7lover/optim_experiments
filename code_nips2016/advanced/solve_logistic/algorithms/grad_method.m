function [xplus, trash] = grad_method(fun,xHist,k,~)

L = fun.L;
gx = fun.fp;
x = xHist(:,k-1);

xplus = x-(1/L)*gx(x);
trash = [];