% solve_logistic_ampe

emptystruct = struct();

%% With LS
warning('off','MATLAB:nearlySingularMatrix')
offrnaParam = [];
offrnaParam.window_size = window_size;
offrnaParam.backtracking = true;
offrnaParam.accelerated = false;
offrnaParam.strong_convex = false;
offrnaParam.online = false;
offrnaParam.lambda = 1e-8;
tic
display('offrna-backt')
[~,offrnabackTolF,offrnabackTolX,~] = do_k_iterations(offrnaParam,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
toffrnaback=toc;
toffrnaback

%% With LS
warning('off','MATLAB:nearlySingularMatrix')
offrnaParam = [];
offrnaParam.window_size = window_size;
offrnaParam.backtracking = true;
offrnaParam.strong_convex = true;
offrnaParam.accelerated = true;
offrnaParam.online = false;
offrnaParam.lambda = 1e-8;
tic
display('nest-offrna-back')
[~,offrnanestbackTolF,offrnanestbackTolX,~] = do_k_iterations(offrnaParam,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
tnestbackoffrna=toc;
tnestbackoffrna

%% 
warning('on','MATLAB:nearlySingularMatrix')
