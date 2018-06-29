% solve_logistic_ampe

emptystruct = struct();

% %% Without LS
% warning('off','MATLAB:nearlySingularMatrix')
% onlinernaParam = [];
% onlinernaParam.window_size = window_size;
% onlinernaParam.backtracking = false;
% onlinernaParam.accelerated = false;
% onlinernaParam.strong_convex = false;
% tic
% display('online-rna')
% [~,ornaTolF,ornaTolTolX,~] = do_k_iterations(onlinernaParam,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
% torna=toc;
% torna


%% With LS
warning('off','MATLAB:nearlySingularMatrix')
onlinernaParam = [];
onlinernaParam.window_size = window_size;
onlinernaParam.backtracking = true;
onlinernaParam.accelerated = false;
onlinernaParam.strong_convex = false;
onlinernaParam.lambda = 1e-18;
tic
display('online-rna-backt')
[~,ornabackTolF,ornabackTolX,~] = do_k_iterations(onlinernaParam,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
tornaback=toc;
tornaback


% %% Without LS
% warning('off','MATLAB:nearlySingularMatrix')
% onlinernaParam = [];
% onlinernaParam.window_size = window_size;
% onlinernaParam.backtracking = false;
% onlinernaParam.strong_convex = true;
% onlinernaParam.accelerated = true;
% tic
% display('nest-rna')
% [~,rnanestTolF,rnanestTolX,~] = do_k_iterations(onlinernaParam,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
% tnestrna=toc;
% tnestrna
% 
% 
% %% With LS
% warning('off','MATLAB:nearlySingularMatrix')
% onlinernaParam = [];
% onlinernaParam.window_size = window_size;
% onlinernaParam.backtracking = true;
% onlinernaParam.strong_convex = true;
% onlinernaParam.accelerated = true;
% tic
% display('nest-rna-back')
% [~,rnanestbackTolF,rnanestbackTolX,~] = do_k_iterations(onlinernaParam,@rna_first_order,nIte,finfo,emptystruct,nIterTol);
% tnestbackrna=toc;
% tnestbackrna

%% 
warning('on','MATLAB:nearlySingularMatrix')
