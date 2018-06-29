function info = rna_sgdlib(problem,param)


k = param.k;
x0 = param.w_init;
n = problem.samples();

k_flat = 0;
k_multiplier = 0;

if isfield(param,'doAdaptiveLambda') % in fact, only 1 data query because we can compute f_i(x) at several datapoints
    if(param.doAdaptiveLambda)
        k_flat = k_flat + 1;
    else
        k_flat = k_flat + 0;
    end
else
    k_flat = k_flat + 1; % default: activated
end


if isfield(param,'forceDecrease')
    if(param.forceDecrease)
        k_flat = k_flat + 1;
    else
        k_flat = k_flat + 0;
    end
else
    k_flat = k_flat + 1; %default: activated
end

epoch_per_iter = k_multiplier * k + k_flat;
epoch_per_iter_method = param.epoch_per_iter_method;
max_iter = floor(param.max_epoch/(k*epoch_per_iter_method+epoch_per_iter));
w = zeros(problem.dim(),max_iter+1);

grad_calc_count = (0:max_iter)*epoch_per_iter;
time = zeros(max_iter+1,1);

w(:,1) = x0;
grad_calc_count(:,1) = 0;
time(1) = 0;

% fill option for sub-algo
if ~isfield(param,'algooption')
    algooption=struct();
else
    algooption = param.algooption;
end
algooption.w_init = x0;
algooption.max_epoch = k;
algooption.store_w = true;
algooption.tol_optgap = 0;
algooption.verbose = -1;
time_beg_tot = cputime;

for iter=2:max_iter+1
%     time_beg = cputime;
    time_beg = tic;
    [~,info_algo] = param.algo(problem, algooption);
    extrapolation = extrapolate(info_algo,problem,param);
%     time_end = cputime;
%     time(iter) = time(iter-1) + time_end-time_beg;
    delta_time = toc(time_beg);
    time(iter) = time(iter-1) + delta_time;
    
    algooption.w_init = extrapolation;
    w(:,iter) = extrapolation;
    grad_calc_count(:,iter) = grad_calc_count(:,iter-1) + info_algo.grad_calc_count(end) + epoch_per_iter * n;
end
time_end_tot = cputime;
time_tot = time_end_tot-time_beg_tot;

if(time(end) > time_tot)
    time = time_tot * time/time(end); % Rescale
end


info.w = w;
info.grad_calc_count = grad_calc_count;
info.time = time;


    function out = extrapolate(info_algo,problem,param)
        finfo.f = @(x) problem.cost(x);
        finfo.proxoperator.f = @(x) 0;
        out = abstract_extrapolation_adaptive_lambda_ls_stepsize(finfo,info_algo.w,param);
    end


end
