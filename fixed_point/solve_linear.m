function [X,tolf,tolX] = solve_linear(A,xstar,x0,k,algoName,L,mu,addParam)
% minimize (x-x^*)'*A*(x-x^*) using k steps of algoName :
% - 'Gradient' for Gradient algorithm
% - 'Nesterov' for Nesterov's algorithm
% - 'AmpeGrad' for ampe acceleration of gradient


param.A = A;
param.x0 = x0;
param.xstar = xstar;

if(nargin < 6)
    param.L = norm(A);
else
    param.L = L;
end

if(nargin<7)
    param.mu = eigs(A,1,'sm');
else
    param.mu = mu;
end

if(nargin<8)
    addParam = [];
end

finfo = getFunction( 'QuadPreDefined' , param);

if(strcmpi('Gradient',algoName))
    grad_param = [];
    [X, tolf, tolX] = do_k_iterations(grad_param,@grad_method_strong_convex,k,finfo);
    
elseif(strcmpi('Nesterov',algoName))
    nest_param.y = finfo.x0;
    [X, tolf, tolX] = do_k_iterations(nest_param,@nesterov_method_strong_convex,k,finfo);
    
elseif(strcmpi('NesterovConvex',algoName))
    nest_param.lambda = 0;
    nest_param.y = finfo.x0;
    [X, tolf, tolX] = do_k_iterations(nest_param,@nesterov_method,k,finfo);
    
elseif(strcmpi('ampeGradient',algoName))
    ampegrad_param.lambda = addParam.lambdaAmpe;
    ampegrad_param.k = addParam.nIteAmpe;
    [X, tolf, tolX] = do_k_iterations(ampegrad_param,@ampe_grad,k,finfo);
    
elseif(strcmpi('ampePower',algoName))
    ampepower_param.lambda = addParam.lambdaAmpe;
    ampepower_param.k = addParam.nIteAmpe;
    [X, tolf, tolX] = do_k_iterations(ampepower_param,@ampe_powerMethod,k,finfo);
    
elseif(strcmpi('ampeNesterov',algoName))
    ampenest_param.lambda = addParam.lambdaAmpe;
    ampenest_param.k = addParam.nIteAmpe;
    [X, tolf, tolX] = do_k_iterations(ampenest_param,@ampe_nest,k,finfo);
    
elseif(strcmpi('conjgrad',algoName))
    [X, tolf, tolX] = conjgrad(finfo,k);
end