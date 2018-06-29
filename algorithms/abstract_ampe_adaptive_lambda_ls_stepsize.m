function [xout,param,fval] = abstract_ampe_adaptive_lambda_ls_stepsize(finfo,x,ite,param)

finfo.x0 = x;

tolfun = NaN; % Avoid computation of tolf

if ~isfield(param,'optialgoparam')
    param.optialgoparam = [];
end

if ~isfield(param,'nPointsToUse')
    param.nPointsToUse = param.k;
end

[algo_x,~,~,optialgoparam] = do_k_iterations(param.optialgoparam,param.optialgo,param.k,finfo,tolfun,param.nPointsToUse);

warning('off','MATLAB:nearlySingularMatrix') % off warnings
warning('off','MATLAB:singularMatrix')
xout = abstract_extrapolation_adaptive_lambda_ls_stepsize(finfo,algo_x,param);
warning('on','MATLAB:singularMatrix')
warning('on','MATLAB:nearlySingularMatrix')


if(any(isnan(xout)))
%     xout = algo_x(:,end);
end
param.optialgoparam = optialgoparam;

if(isfield(param,'useprox') && param.useprox == true)
    if ( isfield(param,'stepsizeprox') )
        xout = finfo.proxoperator.prox(xout,norm(xout-x)*param.stepsizeprox);
    else
        xout = finfo.proxoperator.prox(xout,0);
    end
end
fval = nan;
