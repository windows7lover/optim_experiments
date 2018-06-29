% solve_logistic
%
% This file solves all problems presented in the paper
% https://papers.nips.cc/paper/6267-regularized-nonlinear-acceleration.pdf
%
% Do not forget to add the whole "advanced" foler in your matlab path, i.e. 
% right-click on "advanced" and select:
%   "Add to path > Selected folder and subfolders"
%
% Note: The results may slightly differ from the paper, since some
% optimizations have been done in the code.

%% Add path automatically

addpath(genpath('algorithms/'))
addpath(genpath('data_logistic/'))


%%

clear all
close all
clc

%% Algo parameters

rmpe_lambda0 = 1e-1;
rmpe_lambdamin = 1e-26;
rmpe_k = 5;
% rmpe_optialgo = @grad_method; % Algorithm we want to accelerate
rmpe_optialgo = @grad_method_strong_convex; % Usually it works better

%% Dataset used in the paper

% dataset = 'madelon';
% reg = 1e7; nIte = 1e2;    % Recommended nIte:1e2
% reg = 1e2; nIte = 1e5;    % Recommended nIte:1e5
% reg = 1e-3; nIte = 1e4;   % Recommended nIte:1e4


dataset = 'sonar';
reg = 1e-1; nIte = 1e3;    % Recommended nIte:1e3
% reg = 1e-6; nIte = 1e5;    % Recommended nIte:1e5


% dataset = 'sido0';
% reg = 1e2; nIte = 1e3;    % Recommended nIte:1e3

% Note: for dataset sido0, the result will differ from the paper. Due to a
% bug in the previous experiments, nestBack was too slow.
% The other results are still unchanged.

[ paramFunction.X,paramFunction.y,nFeatures,nPoints,xstar,fstar,paramFunction.lambda ] = load_data( dataset, reg );

%%

paramFunction.x0 = zeros(nFeatures+1,1);
paramFunction.approxL = false;

finfo = getFunction('logistic' , paramFunction);

finfo.xstar = xstar;
finfo.fstar = fstar;

%%

display('Gradient')
tic
[~,tolFgrad,tolXgrad,~] = do_k_iterations([],@grad_method_strong_convex,nIte,finfo);
tgrad = toc;
tgrad

display('Nesterov')
tic
[~,tolFnest,tolXnest,~] = do_k_iterations([],@nesterov_method_strong_convex,nIte,finfo);
tnest=toc;
tnest

display('NesterovBacktracking')
tic
[~,tolFnestback,tolXnestback,~] = do_k_iterations([],@nesterov_method_strong_convex_backtracking,nIte,finfo);
tnextback = toc;
tnextback 

%%

warning('off','MATLAB:nearlySingularMatrix')

RmpeParam.doAdaptiveLambda = true; % Determine if lambda should change over time
RmpeParam.lambda = rmpe_lambda0;
RmpeParam.lambdamin = rmpe_lambdamin; % These two value determine the range of the grid search [lambdamin,lambda0]
RmpeParam.k = rmpe_k; % Number of iteration before extrapolation then restart
RmpeParam.forceDecrease = true; % Optionnal, check if the extrapolated value is smaller than the last iterate of gradient method
RmpeParam.optialgo=rmpe_optialgo; % Algorithm we want to accelerate, here we used by default grad_method_strong_convex
RmpeParam.rmpealgo=@rmpe; % Extrapolation algorithm

nIteRmpe = floor(nIte/RmpeParam.k);

% Without LS
RmpeParam.doLineSearch = false;
tic
display('RMPE')
[~,rmpeGradTolF5,rmpeGradTolX5,~] = do_k_iterations(RmpeParam,@rmpe_adaptive_lambda_ls_stepsize,nIteRmpe,finfo);
trmpe5=toc;
trmpe5


% With LS
RmpeParam.doLineSearch = true;
tic
display('RMPE-LS')
[~,rmpeGradTolF5_ls,rmpeGradTolX5_ls,~] = do_k_iterations(RmpeParam,@rmpe_adaptive_lambda_ls_stepsize,nIteRmpe,finfo);
trmpe5_ls = toc;
trmpe5_ls

warning('on','MATLAB:nearlySingularMatrix')

%% Plot

plot_logistic

