function prox_operator = get_prox(name,param)

if strcmpi(name,'noprox')
    prox_operator.f = @(x) 0;
    prox_operator.fp = @(x) zeros(size(x));
    prox_operator.prox = @(x,h) x;
    prox_operator.mu = 0;
end

if(strcmpi(name,'elasticnet'))
    regl1 = param.regl1;
    regl2 = param.regl2;
    prox_operator.f = @(x) regl1*norm(x,1) + (regl2/2)*norm(x,2)^2;
    prox_operator.fp = @(x) regl1*sign(x) + regl2*x;
    prox_operator.prox = @(x,h) prox_elasticnet(x,h,regl1,regl2);
    prox_operator.mu = regl2;
end