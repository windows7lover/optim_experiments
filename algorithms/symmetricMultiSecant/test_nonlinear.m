% test_nl

N = 5;
nIter = 75;

%% Dataset used in the paper

lambda = 1e-8; lambda_name = 'lambda_1e-8';
% lambda = 0; lambda_name = 'lambda_0';

% dataset = 'sonar'; approxL = false;
dataset = 'madelon'; approxL = false;
% dataset = 'sido0'; approxL = true;

% load dataset
precent_test = 0;
[ paramFunction.X,paramFunction.y,nFeatures,nPoints,paramFunction.Xtest,paramFunction.ytest,nPointsTest] = load_data(dataset, precent_test);


L_dataset = compute_Lmax_logistic(paramFunction.X);

% reg = L_dataset/(nPoints/100); reg_name = 'well_cond';  max_epoch = 100; % well conditionned
% reg = L_dataset/nPoints; reg_name = 'regular_cond';  max_epoch = 100; % normal conditionned
reg = L_dataset/(100*nPoints); reg_name = 'bad_cond';  max_epoch = 100; % badly conditionned


paramFunction.lambda = reg;
paramFunction.approxL = false;

%%

paramFunction.x0 = zeros(nFeatures+1,1);
finfo = getFunction('logistic' , paramFunction);
finfo.L = L_dataset + reg;

solve_problem = true;
% solve_problem = false;


%%

x0 = finfo.x0;
d = length(x0);

f = @(x)finfo.f(x);
g = @(x) finfo.fp(x)/finfo.L;
accum_nlSymQN = Accumulator(d,N);
accum_SymmHessian = Accumulator(d,N);
accum_rna = Accumulator(d,N);


err_fun = @(x) f(x);
% err_fun = @(x) norm(g(x));

err_grad = zeros(nIter,1);
err_bfgs_SymmHessian = zeros(nIter,1);
err_bfgs_nlSymQN = zeros(nIter,1);
err_rna = zeros(nIter,1);


% Gradient
y=x0;
for i=1:nIter
    x = y-g(y);
    err_grad(i) = err_fun(y) ;
    y=x;
    i
end


% BFGS nlSymQN
y=x0;
for i=1:nIter
    x = y-g(y);
    accum_nlSymQN = accum_nlSymQN.store(y,x-y);
    [qNstep] = nonLinearPrecondQN(accum_nlSymQN,1,-g(y));
    err_bfgs_nlSymQN(i) = err_fun(y) ;
    y = y-qNstep;
    i
end


% % BFGS non precond
% y=x0;
% for i=1:nIter
%     x = y-g(y);
%     accum_SymmHessian = accum_SymmHessian.store(y,x-y);
%     [qNstep] = getSymHessian(accum_nlSymQN,-1,-g(y));
%     err_bfgs_SymmHessian(i) = err_fun(y) ;
%     y = y+qNstep;
%     i
% end

% BFGS
param.minFuncOpt = [];
param.minFuncOpt.Method = 'lbfgs';
param.minFuncOpt.Corr = N;
[solution,lbfgs_f,algo_tolx,algoparam,time_ite,fcount_bfgs] = dokiter_minfunc(param,[],nIter,finfo);


% RNA
y=x0;
for i=1:nIter
    x = y-g(y);
    accum_rna = accum_rna.store(y,x-y);
    err_rna(i) = err_fun(y) ;
    y = rna(accum_rna);
    i
end


fstar = min(err_bfgs_nlSymQN)-1e-15;
fstar = min(min(err_rna)-1e-15,fstar);
fstar = min(min(lbfgs_f)-1e-15,fstar);

semilogy(err_grad-fstar)
hold on
semilogy(err_bfgs_nlSymQN-fstar)
semilogy(err_rna-fstar)
semilogy((1:length(lbfgs_f))*fcount_bfgs/length(lbfgs_f), lbfgs_f-fstar)

legend({'gradient','multisecant nlSymQN','rna','lbfgs'})