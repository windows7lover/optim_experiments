% comparison with nesterov

%% Add path automatically
addpath(genpath('../algorithms/'))
addpath(genpath('../classif'))
clear all;
clc;
% close all;
window_size = 10;

backtracking = true;
% backtracking = false;

%% Dataset used in the paper

% dataset = 'sonar'; approxL = false; nIte = 50;
% dataset = 'madelon'; approxL = false; nIte = 100;
dataset = 'sido0'; approxL = false; nIte = 80;

% load dataset
precent_test = 0;
[ paramFunction.X,paramFunction.y,nFeatures,nPoints,paramFunction.Xtest,paramFunction.ytest,nPointsTest] = load_data(dataset, precent_test);

paramFunction.lambda = 0;
paramFunction.approxL = approxL;
load_function

reg = 1e-3*finfo.L; reg_name = 'well_cond';  % well conditionned
% reg = 1e-6*finfo.L; reg_name = 'regular_cond'; % normal conditionned
% reg = 1e-9*finfo.L; reg_name = 'bad_cond'; % badly conditionned


% reg = 1e-2*finfo.L; reg_name = 'comp_nesterov_m4';  % well conditionned
paramFunction.lambda = reg;

load_function;

solve_problem = true;
% solve_problem = false;

%% Algo parameters
nIterTol = NaN;


%%

finfo.mu = 10*finfo.mu ;

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



%%

algo_determ;
algoCell_vanilla = {nest,nest_sc};
algoCell_online = {gradient,nest,nest_sc};
% algoCell = {nest};

nAlgo = length(algoCell);
legendCell = cell(1,nAlgo);
lambda = 1e-8; lambda_name = 'lambda_1e-8';

for i=1:length(algoCell_vanilla)
    algoCell_vanilla{i}.name = ['\textbf{(V)} ', algoCell_vanilla{i}.name];
    algoCell_vanilla{i}.param.backtracking = backtracking;
    algoCell_vanilla{i}.linestyle = '--';
end

for i=1:length(algoCell_online)
    algo = algoCell_online{i};
    param = algo.param;
    
    param.backtracking = backtracking;
    param.dorna = true;
    param.online = true;
    param.window_size = window_size;
    param.lambda = lambda;
    
    algoCell_online{i}.param = param;
    algoCell_online{i}.name = ['\textbf{(O)} ', algoCell_online{i}.name];
    algoCell_online{i}.linestyle = '-';
end


algoCell_all = {algoCell_vanilla{:}, algoCell_online{:}};
warning off
paramalgo = {}
for i=1:length(algoCell_all)
    i
    algo = algoCell_all{i};
    param = algo.param;
    tic;
    [~,fval,~,paramalgo{i}] = algo.algo(param,nIte,finfo,nIterTol);
    time = toc;
    algoCell_all{i}.fval = fval;
    algoCell_all{i}.time = time;
end
warning on

algoCell_all{1}.name = 'Nest (convex)';
algoCell_all{2}.name = 'Nest (str convex)';
algoCell_all{3}.name = 'Online RNA';
algoCell_all{4}.name = 'ORNA+Nest (convex)';
algoCell_all{5}.name = 'ORNA+Nest (str convex)';
% algoCell_all{1}.name = ''
% algoCell_all{1}.name = ''

for i=1:length(algoCell_all)
    legendCell{i} = algoCell_all{i}.name;
end

%% plots
fs = 16;
lw = 2;

figure
for i=1:length(algoCell_all)
    algo = algoCell_all{i};
    iter = (1:length(algo.fval))-1;
    plot_algo(algo,lw,{});
    hold on
end
legend(legendCell,'fontsize',fs,'box','off','location','sw');
set(gca,'fontsize',fs)
axis([0 80 1e-16, 1])
xlabel('Iteration (Gradient oracle calls)','interpreter','latex')
ylabel('$f(x)-f(x^*)$','interpreter','latex')
set(gca,'fontsize',fs,'color','none')
