
%%% Min square

%%

clear all
close all
clc

addpath(genpath('../algorithms'))
addpath(genpath('../fixed_point'))

%% Create Problem: min ||Ax-b||^2 + lambda ||x||^2

type = 'RandomStrongConvex';

opt.L = 1e6;
opt.mu = 1e-6; % could be zero!
lambda = 0; % could be zero too!
n = 50;
% rankA = round(sqrt(n));
rankA = n;
% A = generateRandMatrix(n,type,opt);
A = rand(rankA,n);
A = opt.L*A/norm(A);
b = ones(rankA,1);

param.A = A;
param.xstar = ones(n,1);
param.b = A*param.xstar;
param.x0 = zeros(n,1);

finfo = getFunction('LeastSquare',param);

% do one gradient step for x0 for beauty of the plot
finfo.L=finfo.L*10;
finfo.x0 = finfo.x0 - (1/finfo.L)*finfo.fp(finfo.x0);

std_noise = 1000000;

% finfo.fp = @(x) finfo.fp(x) + std_noise*randn(n,1); % additive noise
finfo.fp = @(x) finfo.fp(x) + std_noise*rand_posdef(n)*(x-param.xstar); % multiplicative noise
finfo.L = (finfo.L + std_noise);

nIteMax = 1000;

%% RMPE param

rmpe_k = 10;

param_rmpe.doAdaptiveLambda = true; % Determine if lambda should change over time
% param_rmpe.lambda = 0; % Determine if lambda should change over time
param_rmpe.lambda = 1; % Determine if lambda should change over time
param_rmpe.lambdamin = 1e-3; % Determine if lambda should change over time
% param_rmpe.lambdaSVD = true; % Determine lambda in function of the SVD
param_rmpe.k = rmpe_k; % Number of iteration before extrapolation then restart
param_rmpe.rmpealgo=@ampe; % Extrapolation algorithm
param_rmpe.doLineSearch = false;
param_rmpe.forceDecrease = false;


%% Declaration for description of algorithms;
nalgo = 2;

% inputs
param_input_cell = cell(nalgo,1);
algohandle_cell = cell(nalgo,1);
nIteMax_cell = cell(nalgo,1);
finfo_cell = cell(nalgo,1);
tolfun_cell = cell(nalgo,1);
nIterTol_cell = cell(nalgo,1);

% ouputs
x_cell = cell(nalgo,1);
tolf_cell = cell(nalgo,1);
tolx_cell = cell(nalgo,1);
param_cell = cell(nalgo,1);
timeite_cell = cell(nalgo,1);

% Plots info
algoname_cell = cell(nalgo,1);
linestyle_cell = cell(nalgo,1);
color_cell = cell(nalgo,1);
plot_ite_multiplicator_cell = cell(nalgo,1);


%% Description of each algo
i = 0;

% Gradient
i = i+1;
param_input_cell{i} = [];
algohandle_cell{i} = @grad_method;
nIteMax_cell{i} = nIteMax;
finfo_cell{i} = finfo;
tolfun_cell{i} = struct();
nIterTol_cell{i} = [];
algoname_cell{i} = 'Gradient';
linestyle_cell{i} = '-.';
color_cell{i} = [1,0,1]; % magenta
plot_ite_multiplicator_cell{i} = 1;


% % Gradient averaging
% i = i+1;
% param_input_cell{i} = [];
% algohandle_cell{i} = @grad_method_averaging;
% nIteMax_cell{i} = nIteMax;
% finfo_cell{i} = finfo;
% tolfun_cell{i} = struct();
% nIterTol_cell{i} = [];
% algoname_cell{i} = 'Gradient Averaging';
% linestyle_cell{i} = '-.';
% color_cell{i} = [0,0,1]; % blue
% plot_ite_multiplicator_cell{i} = 1;

% RMPE - Gradient
i = i+1;
param_rmpe_grad = param_rmpe;
param_rmpe_grad.optialgo=@grad_method;
param_input_cell{i} = param_rmpe_grad;
algohandle_cell{i} = @abstract_ampe_adaptive_lambda_ls_stepsize;
nIteMax_cell{i} = round(nIteMax/param_rmpe_grad.k);
finfo_cell{i} = finfo;
tolfun_cell{i} = struct();
nIterTol_cell{i} = [];
algoname_cell{i} = 'RMPE-gradient';
linestyle_cell{i} = '-';
color_cell{i} = [0,1,0]; % green
plot_ite_multiplicator_cell{i} = param_rmpe_grad.k;


% % RMPE - Gradient averaging
% i = i+1;
% param_rmpe_grad_ave = param_rmpe;
% param_rmpe_grad_ave.optialgo=@grad_method_averaging;
% param_input_cell{i} = param_rmpe_grad_ave;
% algohandle_cell{i} = @abstract_ampe_adaptive_lambda_ls_stepsize;
% nIteMax_cell{i} = round(nIteMax/param_rmpe_grad_ave.k);
% finfo_cell{i} = finfo;
% tolfun_cell{i} = struct();
% nIterTol_cell{i} = [];
% algoname_cell{i} = 'RMPE-gradient averaging';
% linestyle_cell{i} = '-';
% color_cell{i} = [0,1,1]; % cyan
% plot_ite_multiplicator_cell{i} = param_rmpe_grad_ave.k;


%% Run algo
for i=1:nalgo
    i
    tic
    [x_cell{i},tolf_cell{i},tolx_cell{i},param_cell{i},timeite_cell{i}] = ...
        do_k_iterations(param_input_cell{i},algohandle_cell{i},nIteMax_cell{i},finfo_cell{i},tolfun_cell{i},nIterTol_cell{i});
    toc
end

%% Plot algo
figure
for i=1:nalgo
    itevec = ((1:length(tolf_cell{i}))-1)*plot_ite_multiplicator_cell{i};
    semilogy( itevec,tolf_cell{i},linestyle_cell{i},'color',color_cell{i})
    hold on
end
% 
% algoname_cell = [algoname_cell ; algoname_cell];
% 
% for i=1:nalgo
%     cum_mean_algo = cummean(x_cell{i},2);
%     
%     tolf = zeros(size(cum_mean_algo,2),1);
%     for j=1:length(tolf)
%         tolf(j) = finfo.f(cum_mean_algo(:,j));
%     end
%     
%     algoname_cell{nalgo+i} = [algoname_cell{nalgo+i}, ' - cumulative mean'];
%     
%     itevec = ((1:length(tolf_cell{i}))-1)*plot_ite_multiplicator_cell{i};
%     semilogy( itevec,tolf,linestyle_cell{i},'color',color_cell{i})
%     hold on
% end

legend(algoname_cell);
