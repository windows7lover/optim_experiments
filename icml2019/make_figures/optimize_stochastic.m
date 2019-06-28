%% Add path automatically
% addpath(genpath('../algorithms/'))
% addpath(genpath('../classif'))
% clear all;
% clc;

%% Dataset used in the paper

window_size = 10;
lambda = 1e-8; lambda_name = 'lambda_1e-8';
% lambda = 0; lambda_name = 'lambda_0';

dataset = 'sonar'; approxL = false;
% dataset = 'madelon'; approxL = false;
% dataset = 'sido0'; approxL = true;

% load dataset
precent_test = 0;
[ paramFunction.X,paramFunction.y,nFeatures,nPoints,paramFunction.Xtest,paramFunction.ytest,nPointsTest] = load_data(dataset, precent_test);


L_dataset = compute_Lmax_logistic(paramFunction.X);

reg = L_dataset/(nPoints/100); reg_name = 'well_cond';  max_epoch = 100; % well conditionned
% reg = L_dataset/nPoints; reg_name = 'regular_cond';  max_epoch = 100; % normal conditionned
% reg = L_dataset/(100*nPoints); reg_name = 'bad_cond';  max_epoch = 100; % badly conditionned


paramFunction.lambda = reg;
paramFunction.approxL = false;

%%

load_function;
finfo.L = L_dataset + reg;

solve_problem = true;
% solve_problem = false;

%% Algo parameters


% nIterTol = 10; % how much we record
nIterTol = nan; % if nan, record everything

problem = logistic_regression(paramFunction.X, paramFunction.y',paramFunction.Xtest,paramFunction.ytest', reg);


% norm 2 reg
proxOperatorSaga = @(x,stepsize,lambda) x/(1+stepsize*lambda);
proxOperatorKatyusha = @(x,grad,gamma) (1/(1+reg*gamma))*(x-gamma*grad);  % argmin_z 1/(2*gamma) * norm(z-x) + grad*z + psi(z);


batch_size = 1;
w_init = zeros(nFeatures+1,1);

L_logistic = finfo.L_max_sample + reg;
step_size = 1/(3*L_logistic);

%%

if(solve_problem)
    display('Solving with BFGS')
    algoparam.minFuncOpt.Method = 'lbfgs';
    finfo.f = problem.cost;
    finfo.fp = problem.full_grad;
    [xstar,~,~,~] = dokiter_minfunc(algoparam,[],1000,finfo);
    finfo.xstar = xstar;
    fstar = finfo.f(xstar)-(1/(2*reg))*norm(finfo.fp(xstar))^2;
    finfo.fstar = fstar;
    problem.xstar = xstar;
    problem.fstar = fstar;
else
    warning('Solve problem desactivated')
    xstar = nan*finfo.x0;
    finfo.fstar = 0;
    problem.fstar = 0;
    problem.xstar = 0*finfo.x0;
end
%%

algo_stoch

% Convert to add RNA
algoCell_vanilla = algoCell;
nAlgo = 3;
for i=1:nAlgo
    algoCell_vanilla{i}.name = ['\textbf{(V)} ', algoCell_vanilla{i}.name];
end


algoCell_offline = algoCell;
for i=1:nAlgo
    algo = algoCell_offline{i};
    param = algo.param;
    
    param.rna.offline = true;
    param.rna.online = false;
    
    algoCell_offline{i}.param = param;
    algoCell_offline{i}.name = ['\textbf{(P)} ', algoCell_offline{i}.name];
    algoCell_offline{i}.linestyle = '--';
end


algoCell_online = algoCell;
for i=1:nAlgo
    algo = algoCell_online{i};
    param = algo.param;
    
    param.rna.offline = false;
    param.rna.online = true;
    
    algoCell_online{i}.param = param;
    algoCell_online{i}.name = ['\textbf{(O)} ', algoCell_online{i}.name];
    algoCell_online{i}.linestyle = '-';
end

algoCell_all = {algoCell_vanilla{:}, algoCell_offline{:}, algoCell_online{:}, katyusha_algo};


warning off
for i=1:length(algoCell_all)
    algo = algoCell_all{i};
    param = algo.param;
    tic;
    if(i==10) % hack katyusha
        problem.prox = proxOperatorKatyusha;
    else
        problem.prox = 'l2_norm';
    end
    [iter,fval,~,~] = algo.algo(param,problem,nIterTol);
    time = toc;
    algoCell_all{i}.iter = iter;
    algoCell_all{i}.fval = fval;
    algoCell_all{i}.time = time;
    legendCell{i} = algo.name;
end
warning on

%% plots
figure
fs = 16;
lw = 2;

for i=1:length(algoCell_all)
    algo = algoCell_all{i};
    iter = algo.iter;
    semilogy(iter,algo.fval,'linewidth',lw,'color',algo.color,'linestyle',algo.linestyle);
    hold on
end
legend(legendCell,'interpreter','latex','fontsize',fs);
set(gca,'fontsize',fs)


% clear paramFunction
% save(['thesis_figures/matfiles/stochastic_logistic_',dataset,'_',reg_name,'_',lambda_name])



