function xout = abstract_extrapolation_adaptive_lambda_ls_stepsize(finfo,iterations,param)

ampealgo = param.rmpealgo;
doLineSearch = param.doLineSearch;

k = size(iterations,2)-1;

if(~isfield(param,'lambda'))
    param.lambda = 1e-3;
end
lambda = param.lambda; % we use lambda as "lambda0"


if(~isfield(param,'doAdaptiveLambda'))
    param.doAdaptiveLambda = true;
end


if(~isfield(param,'lambdaSVD'))
    param.lambdaSVD = true;
end

if(~isfield(param,'forceDecrease'))
    param.forceDecrease = true;
end

if(~isfield(param,'useBFGS'))
    param.useBFGS = false;
end


if(~isfield(param,'nLambda'))
    param.nLambda = k;
end

algo_x = iterations;
UU = diff(algo_x,1,2);


if(param.useBFGS)
    y = UU(:,1:end-1);
    s = diff(UU,1,2);
    UU = bfgs_multiply([],s,y,UU);
end

UU = UU'*UU;
UUnorml2 = norm(UU);
if(UUnorml2 >= eps)
    UU_norm = UU/UUnorml2;
elseif(isnan(UUnorml2))
    xout = algo_x(:,end);
    return
else
    UU_norm = UU;
end

if(param.doAdaptiveLambda == true)
    if(param.lambdaSVD == true)
        vec = log(svd(UU_norm));
        vec = (vec(1:end-1)+vec(2:end))/2;
        lambdavec=[0 ; exp(vec)];
    else
        lb = param.lambdamin;
        lambdavec = [0, logspace(log10(lb),log10(lambda),param.nLambda)];
    end
else
    lambdavec = lambda;
end

fvalvec = zeros(size(lambdavec));
funval = @(x) finfo.f(x) + finfo.proxoperator.f(x);


for i=1:length(lambdavec)
    xout = ampealgo(algo_x, UU_norm, lambdavec(i));
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
    xout = ampealgo(algo_x, UU_norm, lambdamin);
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
