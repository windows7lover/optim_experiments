% solve_logistic_vanilla_stoch

%%%%
%%%%  HACK : Removed infos.gnorm
%%%%

emptystruct = struct();

%%
display('Saga')
problem_saga = problem;
problem_saga.prox = 'l2_norm';
[w_saga, info_saga] = sag(problem_saga, sgdlib_opt_saga);
[itersaga,tolFsaga,tolXsaga,time_saga,pred_saga] = compute_tol_iterates(info_saga,problem,nIterTol);
time_saga(end)

%%

display('Sgd')
[w_sgd, info_sgd] = sgd(problem, sgdlib_opt_sgd);
[itersgd,tolFsgd,tolXsgd,time_sgd,pred_sgd] = compute_tol_iterates(info_sgd,problem,nIterTol);
time_sgd(end)

%%

display('Svrg')
sgdlib_opt_svrg.max_epoch = sgdlib_opt.max_epoch/2;
[w_svrg, info_svrg] = svrg(problem, sgdlib_opt_svrg);
[itersvrg,tolFsvrg,tolXsvrg,time_svrg,pred_svrg] = compute_tol_iterates(info_svrg,problem,nIterTol);
time_svrg(end)

%%

display('Katyusha')
sgdlib_opt_katyusha.max_epoch = sgdlib_opt.max_epoch/5; % Katyusha compute 5 gradient per iter
problem_katyusha = problem;
problem_katyusha.prox = proxOperatorKatyusha;
[w_katyusha, info_katyusha] = katyusha(problem_katyusha, sgdlib_opt_katyusha);
[iterkatyusha,tolFkatyusha,tolXkatyusha,time_katyusha,pred_katyusha] = compute_tol_iterates(info_katyusha,problem,nIterTol);
time_katyusha(end)