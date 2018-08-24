%% Add path automatically
addpath(genpath('../algorithms/'))
addpath(genpath('../classif'))
clear all;
clc;
% close all;
window_size = 10;

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

% reg = 1e-3*finfo.L; reg_name = 'well_cond';  % well conditionned
reg = 1e-6*finfo.L; reg_name = 'regular_cond'; % normal conditionned
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


%% Solve

algo_determ;
nAlgo = length(algoCell);
legendCell = cell(1,nAlgo);
lambda = 1e-8; lambda_name = 'lambda_1e-8';
% lambda = 0; lambda_name = 'lambda_0';

% Convert to add RNA
algoCell_vanilla = algoCell;
for i=1:nAlgo
    algoCell_vanilla{i}.name = ['\textbf{(V)} ', algoCell_vanilla{i}.name];
end


algoCell_offline = algoCell;
for i=1:nAlgo
    algo = algoCell_offline{i};
    param = algo.param;
    
    param.dorna = true;
    param.online = false;
    param.window_size = window_size;
    param.lambda = lambda;
    
    algoCell_offline{i}.param = param;
    algoCell_offline{i}.name = ['\textbf{(P)} ', algoCell_offline{i}.name];
    algoCell_offline{i}.linestyle = '--';
end


algoCell_online = algoCell;
for i=1:nAlgo
    algo = algoCell_online{i};
    param = algo.param;
    
    param.backtracking = true;
    param.dorna = true;
    param.online = true;
    param.window_size = window_size;
    param.lambda = lambda;
    
    algoCell_online{i}.param = param;
    algoCell_online{i}.name = ['\textbf{(O)} ', algoCell_online{i}.name];
    algoCell_online{i}.linestyle = '-';
end

algoCell_all = {algoCell_vanilla{:}, algoCell_offline{:}, algoCell_online{:}, bfgs_algo};

% vanilla
warning off
for i=1:length(algoCell_all)
    i
    algo = algoCell_all{i};
    param = algo.param;
    tic;
    [~,fval,~,~] = algo.algo(param,nIte,finfo,nIterTol);
    time = toc;
    algoCell_all{i}.fval = fval;
    algoCell_all{i}.time = time;
    legendCell{i} = algo.name;
end
warning on

%% plots
fs = 16;
lw = 2;

figure
for i=1:length(algoCell_all)
    algo = algoCell_all{i};
    iter = (1:length(algo.fval))-1;
    semilogy(iter,algo.fval,'linewidth',lw,'color',algo.color,'linestyle',algo.linestyle);
    hold on
end
legend(legendCell,'interpreter','latex','fontsize',fs);
set(gca,'fontsize',fs)

% save(['thesis_figures/matfiles/logistic_',dataset,'_',reg_name,'_',lambda_name])


