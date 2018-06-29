function [w, infos] = sag(problem, options)
% Stochastic average descent (SAG) algorithm.
%
% Inputs:
%       problem     function (cost/grad/hess)
%       options     options
% Output:
%       w           solution of w
%       infos       information
%
% References:
%       N. L. Roux, M. Schmidt, and F. R. Bach,
%       "A stochastic gradient method with an exponential convergence rate for finite training sets,"
%       NIPS, 2012.
%
% This file is part of SGDLibrary.
%
% Created by H.Kasai on Feb. 15, 2016
% Modified by H.Kasai on Jan. 12, 2017


% set dimensions and samples
d = problem.dim();
n = problem.samples();

% extract options
if ~isfield(options, 'step_init')
    step_init = 0.1;
else
    step_init = options.step_init;
end
step = step_init;

if ~isfield(options, 'step_alg')
    step_alg = 'fix';
else
    if strcmp(options.step_alg, 'decay')
        step_alg = 'decay';
    elseif strcmp(options.step_alg, 'decay-2')
        step_alg = 'decay-2';
    elseif strcmp(options.step_alg, 'fix')
        step_alg = 'fix';
    else
        step_alg = 'decay';
    end
end

if ~isfield(options, 'lambda')
    lambda = 0.1;
else
    lambda = options.lambda;
end

if ~isfield(options, 'tol_optgap')
    tol_optgap = 1.0e-12;
else
    tol_optgap = options.tol_optgap;
end

if ~isfield(options, 'batch_size')
    batch_size = 10;
else
    batch_size = options.batch_size;
end
num_of_bachces = floor(n / batch_size);

if ~isfield(options, 'max_epoch')
    max_epoch = 100;
else
    max_epoch = options.max_epoch;
end

if ~isfield(options, 'w_init')
    w = randn(d,1);
else
    w = options.w_init;
end

if ~isfield(options, 'sub_mode')
    sub_mode = 'SAG';
else
    sub_mode = options.sub_mode;
end

if ~isfield(options, 'f_opt')
    f_opt = -Inf;
else
    f_opt = options.f_opt;
end

if ~isfield(options, 'verbose')
    verbose = false;
else
    verbose = options.verbose;
end

if ~isfield(options, 'store_w')
    store_w = false;
else
    store_w = options.store_w;
end

% initialize
iter = 0;
epoch = 0;
grad_calc_count = 0;

Xtrain = problem.x();


if(~isfield(options,'rna'))
    options.rna = [];
end

if(~isfield(options.rna,'offline'))
    options.rna.offline = false;
end

if(~isfield(options.rna,'online'))
    options.rna.online = false;
end

if(options.rna.offline || options.rna.online)
    rna_module = AccelerationModule(options.rna.window_size,options.rna.lambda);
    if(options.rna.offline && options.rna.online)
        error('Cannot perform offline AND online acceleration at the same time')
    end
else
    options.rna.offline = false;
    options.rna.online = false;
end




perm_idx = randperm(n); % should be the same over time
isScalar = false;

if( strcmpi(problem.name() , 'linear_regression') && batch_size == 1)
    
    grad_array = zeros(1, num_of_bachces);
    isScalar = true;
    
elseif( strcmpi(problem.name() , 'logistic_regression') && batch_size == 1)
    
    grad_array = zeros(1, num_of_bachces);
    isScalar = true;
else
    
    grad_array = zeros(d, num_of_bachces);
    isScalar = false;
    
end

if(isScalar)
    grad_array = zeros(1, num_of_bachces);
    grad_ave = zeros(d, 1);
    
    for j=1:num_of_bachces
        start_index = (j-1) * batch_size + 1;
        indice_j = perm_idx(start_index:start_index+batch_size-1);
        grad_array(1, j) = problem.scal_grad_no_reg(w, indice_j);
        grad_ave = grad_ave + grad_array(1, j)*Xtrain(:,indice_j)/num_of_bachces;
    end
    grad_calc_count = grad_calc_count + num_of_bachces * batch_size;
else
    % prepare an array of gradients, and a valiable of average gradient
    grad_array = zeros(d, num_of_bachces);
    
    for j=1:num_of_bachces
        start_index = (j-1) * batch_size + 1;
        indice_j = perm_idx(start_index:start_index+batch_size-1);
        grad_array(:, j) = problem.grad_no_reg(w, indice_j);
    end
    grad_ave = mean(grad_array, 2);
    grad_calc_count = grad_calc_count + num_of_bachces * batch_size;
end




% store first infos
clear infos;
infos.iter = epoch;
infos.time = 0;
infos.grad_calc_count = grad_calc_count;
%     f_val = problem.cost(w);
%     optgap = f_val - f_opt;
%     infos.optgap = optgap;
%     infos.gnorm = norm(problem.full_grad(w));
%     infos.cost = f_val;
if store_w
    infos.w = w;
end

% display infos
if verbose > 0
    %         fprintf('%s: Epoch = %03d, cost = %.16e, optgap = %.4e\n', sub_mode, epoch, f_val, optgap);
end

% set start time
start_time = tic();

% main loop
optgap = inf;

isSAG = strcmp(sub_mode, 'SAG');
lambda_prob = problem.lambda();

prox_l2 = false;
if(isstr(problem.prox))
    if(strcmpi(problem.prox,'l2_norm'))
        prox_l2 = true;
    end
end


% Main loop

while (optgap > tol_optgap) && (epoch < max_epoch)
    
    batchidx_vec = randi(num_of_bachces,num_of_bachces,1);
    wold = w;
    for j=1:num_of_bachces
        
        batchidx = batchidx_vec(j);
        start_index = (batchidx-1) * batch_size + 1;
        indice_j = perm_idx(start_index:start_index+batch_size-1);
        
        
        if isSAG
            grad = problem.grad(w, indice_j);
            grad_array(:, batchidx) = grad;
            grad_ave = grad_ave + (grad - grad_array(:, batchidx)) / num_of_bachces;
            w = w - step * grad_ave;
        else % SAGA
            
            
            if(isScalar)
                grad = problem.scal_grad_no_reg(w, indice_j);
                
                diff_grad = (grad - grad_array(1, batchidx))*Xtrain(:,indice_j);
                
                grad_array(1, batchidx) = grad;
                
            else
                grad = problem.grad_no_reg(w, indice_j);
                diff_grad = grad - grad_array(:, batchidx);
                grad_array(:, batchidx) = grad;
                
            end
            
            if(prox_l2)
                w = (w - step * ( diff_grad + grad_ave))/(1+step*lambda_prob);
            else
                w = (w - step * ( diff_grad + grad_ave));
                w = problem.prox(w,step,lambda_prob);
            end
            
            grad_ave = grad_ave + diff_grad/num_of_bachces ;
            
            
        end
        
        iter = iter + 1;
    end
    
    
    % measure elapsed time
    elapsed_time = toc(start_time);
    
    % count gradient evaluations
    grad_calc_count = grad_calc_count + num_of_bachces * batch_size;
    % update epoch
    epoch = epoch + 1;
    % calculate optgap
    %         f_val = problem.cost(w);
    %         optgap = f_val - f_opt;
    % calculate norm of full gradient
    % gnorm = norm(problem.full_grad(w));
    
    % store infos
    infos.iter = [infos.iter epoch];
    infos.time = [infos.time elapsed_time];
    infos.grad_calc_count = [infos.grad_calc_count grad_calc_count];
    %         infos.optgap = [infos.optgap optgap];
    %         infos.cost = [infos.cost f_val];
    %         infos.gnorm = [infos.gnorm gnorm];
    
    
    if(options.rna.offline || options.rna.online)
        rna_module = rna_module.store(wold,w,w-wold);
        w_rna = rna_module.accelerate();
    end
        
    if(options.rna.online)
        w = w_rna;
    end
    
    if store_w
        if(options.rna.offline)
            infos.w = [infos.w w_rna];
        else
            infos.w = [infos.w w];
        end
    end
    
    % display infos
    if verbose > 0
        %             fprintf('%s: Epoch = %03d, cost = %.16e, optgap = %.4e\n', sub_mode, epoch, f_val, optgap);
    end
end

if optgap < tol_optgap
    %         fprintf('Optimality gap tolerance reached: tol_optgap = %g\n', tol_optgap);
elseif epoch == max_epoch
    %         fprintf('Max epoch reached: max_epochr = %g\n', max_epoch);
end


if(nargout == 1)
    w = infos;
end

end
