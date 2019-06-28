% launch all experiments

close all
clear all
clc

window_size = 10;
lambda = 1e-8; lambda_name = 'lambda_1e-8';
% lambda = 0; lambda_name = '0';

dataset = 'sonar'; approxL = false;
% dataset = 'madelon'; approxL = false;
% dataset = 'sido0'; approxL = true;

% load dataset
precent_test = 0;
[ paramFunction.X,paramFunction.y,nFeatures,nPoints,paramFunction.Xtest,paramFunction.ytest,nPointsTest] = load_data(dataset, precent_test);


L_dataset = compute_Lmax_logistic(paramFunction.X);
reg = L_dataset/(nPoints/100); reg_name = 'well_cond';  max_epoch = 100; % well conditionned

optimize_stochastic

reg = L_dataset/nPoints; reg_name = 'regular_cond';  max_epoch = 100; % normal conditionned
optimize_stochastic

reg = L_dataset/(100*nPoints); reg_name = 'bad_cond';  max_epoch = 100; % badly conditionned
optimize_stochastic