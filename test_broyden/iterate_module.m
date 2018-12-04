function [y,module,meter] = iterate_module(y,g,nIter,module,finfo,online)

meter = OptimMeter(nIter,finfo);
for i=1:nIter
%     fpold = finfo.fp(y);
    x = g(y);
    module = module.store(y,x);
    extrapolation = module.accelerate();
    meter = meter.store(extrapolation);
    if(online)
        y = extrapolation;
%         fpnew = finfo.fp(y);
%         fpold'*fpnew
    else
        y = x;
    end
end