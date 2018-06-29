function [xplus, trash] = grad_method_strong_convex(fun,xHist,k,~)

L = fun.L;
mu = fun.mu;
gx = fun.fp;
x = xHist(:,k-1);

xplus = x-(2/(mu+L))*gx(x);
trash = [];

