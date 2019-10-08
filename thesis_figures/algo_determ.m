% solve

emptystruct = struct();

color_grad = 1;
color_nest = 2;
color_nest_sc = 3;
color_lbfgs = 4;

nAlgo = 4;
colors = linspecer(nAlgo);


if(exist('ls','var'))
    ls = backtracking
else
    ls = true;
    ls = false;
end

% Define Gradient Method
param.dorna = false;
param.online = false;
param.window_size = 0;
param.backtracking = ls;
param.accelerated = false;
param.strong_convex = false;

gradient.param = param;
gradient.algo = @(param,nIte,finfo,nIterTol) do_k_iterations(param,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
gradient.name = 'Gradient';
gradient.color = colors(color_grad,:);
gradient.linestyle = ':';


% Define Nesterov Method (not strongly convex)
param.dorna = false;
param.online = false;
param.window_size = 0;
param.backtracking = ls;
param.accelerated = true;
param.strong_convex = false;

nest.param = param;
nest.algo = @(param,nIte,finfo,nIterTol) do_k_iterations(param,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
nest.name = 'Nesterov ($1/k$)';
nest.color = colors(color_nest,:);
nest.linestyle = ':';


% Define Nesterov Method (strongly convex)
param.dorna = false;
param.online = false;
param.window_size = 0;
param.backtracking = ls;
param.accelerated = true;
param.strong_convex = true;

nest_sc.param = param;
nest_sc.algo = @(param,nIte,finfo,nIterTol) do_k_iterations(param,@rna_first_order,10*nIte,finfo,emptystruct,nIterTol);
nest_sc.name = 'Nesterov ($1-\sqrt{\kappa}$) ';
nest_sc.color = colors(color_nest_sc,:);
nest_sc.linestyle = ':';

algoCell = {gradient,nest,nest_sc};



% Define lBFGS
param.dorna = false;
param.online = false;
param.window_size = 0;
param.backtracking = ls;
param.accelerated = true;
param.strong_convex = true;

nest_sc.param = param;
nest_sc.algo = @(param,nIte,finfo,nIterTol) do_k_iterations(param,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
nest_sc.name = 'Nesterov ($1-\sqrt{\kappa}$) ';
nest_sc.color = colors(color_nest_sc,:);
nest_sc.linestyle = ':';


param.minFuncOpt = [];
param.minFuncOpt.Method = 'lbfgs';
param.minFuncOpt.Corr = window_size;
bfgs_algo.param = param;
bfgs_algo.algo = @(param,nIte,finfo,nIterTol) dokiter_minfunc(param,[],nIte,finfo);
bfgs_algo.name = 'l-BFGS';
bfgs_algo.color = colors(color_lbfgs,:);
bfgs_algo.linestyle = '-';
