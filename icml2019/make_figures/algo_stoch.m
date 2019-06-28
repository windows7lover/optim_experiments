% algo_stoch% solve

emptystruct = struct();

color_sgd = 1;
color_svrg = 2;
color_saga = 3;
color_katyusha = 4;

nAlgo = 4;
colors = linspecer(nAlgo);


%%% Options for vanilla algo
L_logistic = finfo.L_max_sample + reg;
stepsize_saga = 1/(3*L_logistic);

sgdlib_opt.max_epoch = max_epoch;    
sgdlib_opt.store_w = true;
sgdlib_opt.tol_optgap = 0;
sgdlib_opt.batch_size = 1;
sgdlib_opt.w_init = w_init;
sgdlib_opt.step_init = stepsize_saga;
sgdlib_opt.rna.window_size = window_size;
sgdlib_opt.rna.lambda = lambda;


% 
% problem_saga = problem;
% problem_saga.prox = 'l2_norm';
% 
% sgdlib_opt_katyusha.L_max_sample = finfo.L_max_sample;
% sgdlib_opt_katyusha.mu = finfo.mu;
% 
% sgdlib_opt_saga.sub_mode = 'saga';



% Define SGD
sgd_algo.param = sgdlib_opt;
sgd_algo.algo = @(param,prob,nIterTol) (compute_tol_iterates(sgd(prob, param),prob ,nIterTol));
sgd_algo.name = 'SGD';
sgd_algo.color = colors(color_sgd,:);
sgd_algo.linestyle = ':';


% Define SVRG

svrg_algo.param = sgdlib_opt;
svrg_algo.param.max_epoch = max_epoch/2;
svrg_algo.algo = @(param,prob,nIterTol) (compute_tol_iterates(svrg(prob, param),prob ,nIterTol));
svrg_algo.name = 'SVRG';
svrg_algo.color = colors(color_svrg,:);
svrg_algo.linestyle = ':';


% Define Saga

saga_algo.param = sgdlib_opt;
saga_algo.param.sub_mode = 'saga';
saga_algo.algo = @(param,prob,nIterTol) (compute_tol_iterates(sag(prob, param),prob ,nIterTol));
saga_algo.name = 'Saga';
saga_algo.color = colors(color_saga,:);
saga_algo.linestyle = ':';


algoCell = {sgd_algo,svrg_algo,saga_algo};

katyusha_algo.param = sgdlib_opt;
katyusha_algo.param.max_epoch = max_epoch/5;
katyusha_algo.param.L_max_sample = finfo.L_max_sample;
katyusha_algo.param.mu = finfo.mu;
katyusha_algo.algo = @(param,prob,nIterTol) (compute_tol_iterates(katyusha(prob, param),prob ,nIterTol));
katyusha_algo.name = 'Katyusha';
katyusha_algo.color = colors(color_katyusha,:);
katyusha_algo.linestyle = '-';
