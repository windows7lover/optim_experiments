% test

%% Add path automatically
addpath(genpath('../algorithms/'))
addpath(genpath('../classif'))
clear all;
clc;

%% Dataset used in the paper

dataset = 'sonar'; approxL = false; nIte = 500;
% dataset = 'madelon'; approxL = false; nIte = 500;
% dataset = 'sido0'; approxL = true; nIte = 500;

% load dataset
precent_test = 0;
[ paramFunction.X,paramFunction.y,nFeatures,nPoints,paramFunction.Xtest,paramFunction.ytest,nPointsTest] = load_data(dataset, precent_test);

paramFunction.lambda = 0;
paramFunction.approxL = approxL;
load_function

reg = 1e-3*finfo.L; reg_name = 'well_cond';  % well conditionned
% reg = 1e-6*finfo.L; reg_name = 'regular_cond'; % normal conditionned
% reg = 1e-9*finfo.L; reg_name = 'bad_cond'; % badly conditionned
paramFunction.lambda = reg;

load_function;

solve_problem = true;
% solve_problem = false;

%% Algo parameters
nIterTol = NaN;


%%

if(solve_problem)
    display('Solving with BFGS')
    algoparam.minFuncOpt.Method = 'lbfgs';
%     if(~determ)
%         finfo.f = problem.cost;
%         finfo.fp = problem.full_grad;
%     end
    [xstar,~,~,~] = dokiter_minfunc(algoparam,[],1000,finfo);
    finfo.xstar = xstar;
    fstar = finfo.f(xstar)-(1/(2*reg))*norm(finfo.fp(xstar))^2;
    finfo.fstar = fstar;
    problem.xstar = xstar;
    problem.fstar = fstar;
end

nIter = 1000;
fseq_adaptive = adapaccel( finfo.x0, finfo, nIter );
fseq_nesterov = nesterov( finfo.x0, finfo, nIter );

figure
semilogy(0:nIter, fseq_adaptive, 'linewidth',2)
hold on
semilogy(0:nIter, fseq_nesterov, 'linewidth',2)
kappa = finfo.mu/finfo.L;
semilogy(0:nIter, fseq_adaptive(1) * (1-sqrt(kappa)).^(0:nIter) )
semilogy(0:nIter, fseq_adaptive(1) * 1./(0:nIter).^2 )

