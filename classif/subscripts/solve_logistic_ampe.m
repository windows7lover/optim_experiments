% solve_logistic_ampe

emptystruct = struct();

warning('off','MATLAB:nearlySingularMatrix')

RmpeParam.doAdaptiveLambda = true; % Determine if lambda should change over time
RmpeParam.lambdaSVD = true; % Determine lambda in function of the SVD
RmpeParam.k = rmpe_k; % Number of iteration before extrapolation then restart
% RmpeParam.nPointsToUse = nPointsToUse; % Number of iteration before extrapolation then restart
RmpeParam.optialgo=rmpe_optialgo; % Algorithm we want to accelerate, here we used by default grad_method_strong_convex
RmpeParam.rmpealgo=@ampe; % Extrapolation algorithm

nIteRmpe = floor(nIte/(RmpeParam.k*RmpeParam.optialgoparam.batchsize));

%% Without LS
RmpeParam.doLineSearch = false;
tic
display('RMPE')
[~,rmpeGradTolF,rmpeGradTolX,~] = do_k_iterations(RmpeParam,@abstract_ampe_adaptive_lambda_ls_stepsize,nIteRmpe,finfo,emptystruct,nIterTol);
trmpe=toc;
trmpe


%% With LS
RmpeParam.doLineSearch = true;
tic
display('RMPE-LS')
[~,rmpeGradTolF_ls,rmpeGradTolX_ls,~] = do_k_iterations(RmpeParam,@abstract_ampe_adaptive_lambda_ls_stepsize,nIteRmpe,finfo,emptystruct,nIterTol);
trmpe_ls = toc;
trmpe_ls

warning('on','MATLAB:nearlySingularMatrix')