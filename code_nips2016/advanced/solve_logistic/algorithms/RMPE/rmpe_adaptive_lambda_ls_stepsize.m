function [xout,param] = rmpe_adaptive_lambda_ls_stepsize(finfo,xHist,ite,param)

x = xHist(:,ite-1);

finfo.x0 = x;

tolfun = NaN; % Avoid computation of tolf

if ~isfield(param,'optialgoparam')
    param.optialgoparam = [];
end

[algo_x,~,~,optialgoparam] = do_k_iterations(param.optialgoparam,param.optialgo,param.k,finfo,tolfun);
xout = extrapolation_adaptive_lambda_ls_stepsize(finfo,algo_x,param);

param.optialgoparam = optialgoparam;
