%% Add path automatically
addpath(genpath('../algorithms/'))
addpath(genpath('../classif'))
clear all;
clc;
% close all;

%% Dataset used in the paper

dataset = 'sonar';
% dataset = 'madelon';
% dataset = 'sido0';

% load dataset
precent_test = 0;
[ paramFunction.X,paramFunction.y,nFeatures,nPoints,paramFunction.Xtest,paramFunction.ytest,nPointsTest] = load_data(dataset, precent_test);

paramFunction.lambda = 0;


% paramFunction.x0 = zeros(nFeatures+1,1);
% paramFunction.x0 = ones(nFeatures+1,1);
paramFunction.x0 = rand(nFeatures+1,1);


getf = @(x) getFunction('logistic' , x); solve_problem = true;
% getf = @(x) getFunction('leastsquare2' , x); solve_problem = false;


finfo = getf(paramFunction);

reg = 1e-3*finfo.L; reg_name = 'well_cond';  % well conditionned
% reg = 1e-6*finfo.L; reg_name = 'regular_cond'; % normal conditionned
% reg = 1e-9*finfo.L; reg_name = 'bad_cond'; % badly conditionned
paramFunction.lambda = reg;

finfo = getf(paramFunction);
if(~isfield(finfo,'lsfun'))
    finfo.lsfun = nan;
end
if(~isfield(finfo,'lsnormgrad'))
    finfo.lsnormgrad = nan;
end


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

warning off

g = @(x) x-finfo.fp(x)/finfo.L;
nIter = 1000;
N = 50;
beta = 1;
% do_online = false;
do_online = true;

lshandle = finfo.lsfun;
% lshandle = finfo.lsnormgrad;


% module_name = {'none','good_anderson','bad_anderson','good_broyden','bad_broyden','bfgs','dfp','conjgrad','gmres'};
module_name = {'none','good_anderson','bad_anderson','good_broyden','bad_broyden','bfgs','dfp'};
module_cell = cell(length(module_name),1);

param.beta = beta;
param.H0 = beta;
param.M = 1;
param.lambda = 1e-6;
% param.linesearch = @finfo.lsfun;

for i=1:length(module_name)
    module_cell{i} = AccelerationModule2(N,module_name{i},param);
end


x0 = finfo.x0;


meter_cell = cell(length(module_cell),1);
for i=1:length(module_cell)
    [~,module_cell{i},meter_cell{i}] = iterate_module(x0,g,nIter,module_cell{i},finfo,do_online);
end

%%
lw = 2;
color = linspecer(length(meter_cell));
for i=1:length(meter_cell)
    meter_cell{i} = meter_cell{i}.plotParam(color(i,:),lw,'-');
end


figure
subplot(2,3,1)
plotName = 'valf';
for i=1:length(meter_cell)
    meter_cell{i}.plot(plotName);
    hold on
end
title('$f(x)-f(x^*)$','interpreter','latex')
axis square


subplot(2,3,2)
plotName = 'normgrad';
for i=1:length(meter_cell)
    meter_cell{i}.plot(plotName);
    hold on
end
title('$\|\nabla f(x)\|$','interpreter','latex')
axis square


subplot(2,3,3)
plotName = 'normx';
for i=1:length(meter_cell)
    meter_cell{i}.plot(plotName);
    hold on
end
title('$\|x-x^*\|$','interpreter','latex')
axis square


subplot(2,3,5)
plotName = 'normx';
for i=1:length(meter_cell)
    meter_cell{i}.plot(plotName);
    hold on
end


legend_cell = cell(length(module_cell),1);
for i=1:length(module_cell)
    legend_cell{i} = module_name{i};
end
legend(legend_cell,'interpreter','none')
title('Legend','interpreter','latex')
axis off

