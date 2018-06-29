%% Add path automatically
addpath(genpath('../algorithms/'))
addpath(genpath('../classif'))
clear all; 
clc;
% close all;

%%
% determ = true;
determ = false;

%% Dataset used in the paper

% dataset = 'sonar'; approxL = false;
dataset = 'madelon'; approxL = false;
% dataset = 'sido0'; approxL = true;

% load dataset
precent_test = 0;
[ paramFunction.X,paramFunction.y,nFeatures,nPoints,paramFunction.Xtest,paramFunction.ytest,nPointsTest] = load_data(dataset, precent_test);


% % generate random dataset
% dataset = 'random';
% nPoints_random = 2000;
% nFeatures_random = 1500;
% precent_test = 0;
% [ paramFunction.X,paramFunction.y,nFeatures,nPoints,paramFunction.Xtest,paramFunction.ytest,nPointsTest] = load_data(dataset, precent_test,nPoints_random,nFeatures_random);


if(~determ)
    L_dataset = compute_Lmax_logistic(paramFunction.X);

%     reg = L_dataset/(nPoints/100); reg_name = 'well_cond';  max_epoch = 150; % well conditionned
%     reg = L_dataset/nPoints; reg_name = 'regular_cond';  max_epoch = 100; % normal conditionned
    reg = L_dataset/(100*nPoints); reg_name = 'bad_cond';  max_epoch = 100; % badly conditionned
else
    paramFunction.lambda = 0;
    paramFunction.approxL = approxL;
    load_function
    
%     reg = 1e-3*finfo.L; reg_name = 'well_cond'; nIte = 50; % well conditionned
%     reg = 1e-6*finfo.L; reg_name = 'regular_cond'; nIte = 50; % normal conditionned
    reg = 1e-9*finfo.L; reg_name = 'bad_cond'; nIte = 50; % badly conditionned
end

paramFunction.lambda = reg;

if(~determ)
   paramFunction.approxL = true;
else
    paramFunction.approxL = false;
end
%%

load_function;

if(~determ)
    finfo.L = L_dataset + reg;
end

solve_problem = true;
% solve_problem = false;

%% Algo parameters


% nIterTol = 10; % how much we record
nIterTol = nan; % if nan, record everything
    
if(determ)
    window_size = 10;
else
    
    
    %%%
    
    % BIG PROBLEM WITH STEP SIZE USED IN THE TOOLBOX!
    % done
    
    % BIG PROBLEM IN COMPUTATION OF FGRAD IN THE TOOLBOX
    % done
    
    % BIG PROBLEM WITH KATYUSHA
    % the regularizer is splitted, we have to remove it and compute prox in katyusha.
    % done
    
    %%%
    
    problem_type = 'leastsquare';
%     problem_type = 'logistic';
       
    if(strcmpi(problem_type,'logistic'))
        problem = logistic_regression(paramFunction.X, paramFunction.y',paramFunction.Xtest,paramFunction.ytest', reg); 
    end
    
    if(strcmpi(problem_type,'leastsquare'))
        problem = linear_regression(paramFunction.X, paramFunction.y',paramFunction.Xtest,paramFunction.ytest', reg);
    end
    
    % norm 2 reg
    proxOperatorSaga = @(x,stepsize,lambda) x/(1+stepsize*lambda);
    proxOperatorKatyusha = @(x,grad,gamma) (1/(1+reg*gamma))*(x-gamma*grad);  % argmin_z 1/(2*gamma) * norm(z-x) + grad*z + psi(z);
    
    w_init = zeros(nFeatures+1,1);
    
end

%%

if(solve_problem)
    display('Solving with BFGS')
    algoparam.minFuncOpt.Method = 'lbfgs';
    if(~determ)
        finfo.f = problem.cost;
        finfo.fp = problem.full_grad;
    end
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
if(determ)
    solve_logistic_vanilla_deterministic
else
    solve_logistic_vanilla_stoch
end

%%
warning off
if(determ)
%     solve_logistic_ampe
    solve_logistic_offline_rna
    solve_logistic_online_rna
else
    solve_logistic_ampe_stoch
end
warning on

% %% Plot
if(determ)
    plot_logistic
else
    name_fig = [ problem_type , '_' , dataset , '_' , reg_name , '_' , num2str(nPoints) , '_' , num2str(nFeatures) ];
    plot_logistic_stoch
%     plot_logistic_stoch_accuracy
end
