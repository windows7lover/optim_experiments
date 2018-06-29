function xout = extrapolation_adaptive_lambda_ls_stepsize(finfo,iterations,param)

rmpealgo = param.rmpealgo;
doLineSearch = param.doLineSearch;

k = size(iterations,2)-1;
lambda = param.lambda; % we use lambda as "lambda0"

if(~isfield(param,'doAdaptiveLambda'))
    param.doAdaptiveLambda = true;
end

if(~isfield(param,'forceDecrease'))
    param.forceDecrease = true;
end

if(param.doAdaptiveLambda == true)
    lb = param.lambdamin;
    lambdavec = [0, logspace(log10(lb),log10(lambda),k)];
else
    lambdavec = lambda;
end

fvalvec = zeros(size(lambdavec));
funval = @(x) finfo.f(x);
algo_x = iterations;
UU = diff(algo_x,1,2);
UU = UU'*UU;
UU_norm = UU/norm(UU);


for i=1:length(lambdavec)
    xout = rmpealgo(algo_x, UU_norm, lambdavec(i));
    if(param.doAdaptiveLambda)
        fvalvec(i) = funval(xout);
    else
        fvalvec = 0;
    end
end

if( param.forceDecrease && min(fvalvec) > funval(algo_x(:,end)) ) % Force decrease
    xout = algo_x(:,end);
else
    [~, idx_min] = min(fvalvec);
    lambdamin = lambdavec(idx_min);
    xout = rmpealgo(algo_x, UU_norm, lambdamin);
end




if(doLineSearch)
    x = iterations(:,1);
    step = xout-x;
    sizestep = 1;
    fold = funval(x+sizestep*step);
    sizestep = 2*sizestep;

    while(true)
        fnew = funval(x+sizestep*step);
        if(fold>fnew)
            fold = fnew;
            sizestep = 2*sizestep;
        else
            break;
        end
    end
    sizestep = sizestep/2;
    xout = x+sizestep*step;
end
