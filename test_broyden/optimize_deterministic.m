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


% getf = @(x) getFunction('logistic' , x); solve_problem = true;
getf = @(x) getFunction('leastsquare2' , x); solve_problem = false;


finfo = getf(paramFunction);

% reg = 1e-3*finfo.L; reg_name = 'well_cond';  % well conditionned
% reg = 1e-6*finfo.L; reg_name = 'regular_cond'; % normal conditionned
reg = 1e-9*finfo.L; reg_name = 'bad_cond'; % badly conditionned
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
nIter = 100;
N = 10;
C = nan;
beta = -1;
% beta = nan;
gamma = @(N) [zeros(N-1,1) ; 1] ;
% gamma = @(N) ones(N,1)/N ;
normalization = 'inverse';
% normalization = 'none';
% do_online = false;
do_online = true;

lshandle = finfo.lsfun;
% lshandle = finfo.lsnormgrad;

module_cell = {};

i=1;
module_cell{i} = AccelerationModule(1,0,'none',beta,normalization,gamma,C,lshandle,finfo); i=i+1;
module_cell{i} = AccelerationModule(N,0,'good_anderson',beta,normalization,gamma,C,lshandle,finfo); i=i+1;
module_cell{i} = AccelerationModule(N,0,'bad_anderson',beta,normalization,gamma,C,lshandle,finfo); i=i+1;
module_cell{i} = AccelerationModule(N,0,'good_broyden',beta,normalization,gamma,C,lshandle,finfo); i=i+1;
module_cell{i} = AccelerationModule(N,0,'bad_broyden',beta,normalization,gamma,C,lshandle,finfo); i=i+1;
module_cell{i} = AccelerationModule(N,0,'bfgs',beta,normalization,gamma,C,lshandle,finfo); i=i+1;
module_cell{i} = AccelerationModule(N,0,'conjgrad',beta,normalization,gamma,C,lshandle,finfo); i=i+1;
% module_cell{i} = AccelerationModule(N,0,'adagrad',beta,normalization,gamma,C,lshandle,finfo); i=i+1;
% bfgsModule = AccelerationModule(N,0,'bfgs',0,normalization,gamma,,C,lshandle,finfo);


x0 = finfo.x0;


meter_cell = cell(length(module_cell),1);
for i=1:length(module_cell)
    [~,module_cell{i},meter_cell{i}] = iterate_module(x0,g,nIter,module_cell{i},finfo,do_online);
%     module_adagrad = AccelerationModule(nIter,0,'adagrad',beta,normalization,gamma,C,lshandle,finfo);
%     [~,module_cell{i},meter_cell{i}] = iterate_module_adagrad(x0,g,module_adagrad,nIter,module_cell{i},finfo,do_online);
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
    legend_cell{i} = module_cell{i}.accelerationType;
end
legend(legend_cell,'interpreter','none')
title('Legend','interpreter','latex')
axis off

