function [xplus,trash] = proximalgradient(fun,x,k,~)

L = fun.L;
gx = fun.fp;

h = 1/L;
gradstep = x-h*gx(x);
xplus = fun.proxoperator.prox(gradstep,h);

trash = [];