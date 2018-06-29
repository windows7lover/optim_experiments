% solve_logistic_vanilla_deterministic

emptystruct = struct();


% %% Without LS
% warning('off','MATLAB:nearlySingularMatrix')
% param = [];
% param.dorna = false;
% param.window_size = 0;
% param.backtracking = false;
% param.accelerated = false;
% param.strong_convex = false;
% tic
% display('Gradient')
% [~,tolFgrad,tolXgrad,~] = do_k_iterations(param,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
% tgrad=toc;
% tgrad


%% With LS
warning('off','MATLAB:nearlySingularMatrix')
param = [];
param.dorna = false;
param.window_size = 0;
param.backtracking = true;
param.accelerated = false;
param.strong_convex = false;
tic
display('Gradient')
[~,tolFgradback,tolXgradback,~] = do_k_iterations(param,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
tgradback=toc;
tgradback

%%

% display('Nesterov')
% tic
% [~,tolFnest,tolXnest,~] = do_k_iterations([],@nesterov_method_strong_convex,nIte,finfo,emptystruct,nIterTol);
% tnest=toc;
% tnest


warning('off','MATLAB:nearlySingularMatrix')
param = [];
param.dorna = false;
param.window_size = 0;
param.backtracking = true;
param.accelerated = true;
param.strong_convex = true;
tic
display('NesterovBacktracking')
[~,tolFnestback,tolXnestback,~] = do_k_iterations(param,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
tnextback=toc;
tnextback


%%

display('LBFGS')
tic
algoparam.minFuncOpt = [];
% algoparam.LS_init = 3;
algoparam.minFuncOpt.Method = 'lbfgs';
% algoparam.minFuncOpt.Corr = 25;
% For the Wolfe line-search, these interpolation strategies are available ('LS_interp'):
%   - 0 : Step Size Doubling and Bisection
%   - 1 : Cubic interpolation/extrapolation using new function and gradient values (default)
%   - 2 : Mixed quadratic/cubic interpolation/extrapolation
% algoparam.minFuncOpt.LS_interp = 2;

% LS_type = 1;
% LS_interp = 2;

[~,tolFlbfgs,~,algoparamlbfgs] = dokiter_minfunc(algoparam,[],nIte,finfo);
tlbfgs = toc;
tlbfgs


%% 
display('LBFGS low mem')
tic
algoparam.minFuncOpt = [];
% algoparam.LS_init = 3;
algoparam.minFuncOpt.Method = 'lbfgs';
algoparam.minFuncOpt.Corr = window_size;
% For the Wolfe line-search, these interpolation strategies are available ('LS_interp'):
%   - 0 : Step Size Doubling and Bisection
%   - 1 : Cubic interpolation/extrapolation using new function and gradient values (default)
%   - 2 : Mixed quadratic/cubic interpolation/extrapolation
% algoparam.minFuncOpt.LS_interp = 2;

% LS_type = 1;
% LS_interp = 2;

[~,tolFlbfgslowmem,~,algoparamlbfgslowmem] = dokiter_minfunc(algoparam,[],nIte,finfo);
tlbfgslowmem = toc;
tlbfgslowmem