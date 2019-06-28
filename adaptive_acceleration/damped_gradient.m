function [x, L_k, fxk, gxk] = damped_gradient(x0, finfo, L_k)

% L_k = L_k/2; % makes it adaptive
x = x0-finfo.fp(x0)/L_k;
newgrad = finfo.fp(x);

while newgrad'*(x0-x) < (norm(newgrad)^2)/L_k
    % display('wtf')
    L_k = L_k*2;
    x = x0-finfo.fp(x0)/L_k;
    newgrad = finfo.fp(x);
end

fxk = finfo.f(x);
gxk = newgrad;
