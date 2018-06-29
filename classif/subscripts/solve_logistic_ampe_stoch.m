% solve_logistic_ampe

emptystruct = struct();

warning('off','MATLAB:nearlySingularMatrix')

RmpeParam.rmpealgo=@ampe; % Extrapolation algorithm


%% SAGA

display('RMPE SAGA')
tic
RmpeParam_saga = RmpeParam;
RmpeParam_saga.algo = @(x,y) sag(x,y) ;
RmpeParam_saga.algooption = algooption_saga;
RmpeParam_saga.epoch_per_iter_method = 1;
RmpeParam_saga.k = rmpe_k_saga; % Number of iteration before extrapolation then restart
RmpeParam_saga.max_epoch = RmpeParam.max_epoch;

info_rmpe_saga = rna_sgdlib(problem_saga,RmpeParam_saga);
[iterrmpe_saga,tolFrmpe_saga,tolXrmpe_saga,time_rmpe_saga,pred_rmpe_saga] = compute_tol_iterates(info_rmpe_saga,problem,nIterTol);
time_rmpe_saga(end)
toc


%% SGD

display('RMPE SGD')

RmpeParam_sgd = RmpeParam;
RmpeParam_sgd.algo = @(x,y) sgd(x,y) ;
RmpeParam_sgd.algooption = algooption_sgd;
RmpeParam_sgd.k = rmpe_k_sgd; % Number of iteration before extrapolation then restart
RmpeParam_sgd.epoch_per_iter_method = 1;
RmpeParam_sgd.max_epoch = RmpeParam.max_epoch;

info_rmpe_sgd = rna_sgdlib(problem,RmpeParam_sgd);
[iterrmpe_sgd,tolFrmpe_sgd,tolXrmpe_sgd,time_rmpe_sgd,pred_rmpe_sgd] = compute_tol_iterates(info_rmpe_sgd,problem,nIterTol);
time_rmpe_sgd(end)


%% SVRG

display('RMPE SVRG')

RmpeParam_svrg = RmpeParam;
RmpeParam_svrg.algo = @(x,y) svrg(x,y) ;
RmpeParam_svrg.algooption = algooption_svrg;
RmpeParam_svrg.k = rmpe_k_svrg; % Number of iteration before extrapolation then restart
RmpeParam_svrg.epoch_per_iter_method = 2;
RmpeParam_svrg.max_epoch = RmpeParam.max_epoch; 

info_rmpe_svrg = rna_sgdlib(problem,RmpeParam_svrg);
[iterrmpe_svrg,tolFrmpe_svrg,tolXrmpe_svrg,time_rmpe_svrg,pred_rmpe_svrg] = compute_tol_iterates(info_rmpe_svrg,problem,nIterTol);
time_rmpe_svrg(end)


%% Katyusha

display('RMPE Katyusha')

RmpeParam_katyusha = RmpeParam;
RmpeParam_katyusha.algo = @(x,y) katyusha(x,y) ;
RmpeParam_katyusha.algooption = algooption_katyusha;
RmpeParam_katyusha.k = rmpe_k_katyusha; % Number of iteration before extrapolation then restart
RmpeParam_katyusha.epoch_per_iter_method = 5; % katyusha compute 5 gradient per iter
RmpeParam_katyusha.max_epoch = RmpeParam.max_epoch; 
problem_katyusha = problem;
problem_katyusha.prox = @(x,grad,gamma) (1/(1+reg*gamma))*(x-gamma*grad);  % argmin_z 1/(2*gamma) * norm(z-x) + grad*z + psi(z);

info_rmpe_katyusha = rna_sgdlib(problem_katyusha,RmpeParam_katyusha);
[iterrmpe_katyusha,tolFrmpe_katyusha,tolXrmpe_katyusha,time_rmpe_katyusha,pred_rmpe_katyusha] = compute_tol_iterates(info_rmpe_katyusha,problem,nIterTol);
time_rmpe_katyusha(end)


%%

warning('on','MATLAB:nearlySingularMatrix')
