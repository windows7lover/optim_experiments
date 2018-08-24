%% Add path automatically
addpath(genpath('../algorithms/'))
addpath(genpath('../classif'))
clear all;
clc;
% close all;

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


grad_step = @(x) x-(1/finfo.L)*finfo.fp(x);

batchSize = 5;
lambda = 1e-6;
rna = AccelerationModule(batchSize,lambda,'rna');
broyden = AccelerationModule(batchSize,lambda,'block_broyden');

x = randn(size(finfo.x0));

for i=1:batchSize
    xplus = grad_step(x);
    rna = rna.store(x,xplus,xplus-x);
    broyden = broyden.store(x,xplus,xplus-x);
    x=xplus;
end

norm(rna.accelerate()-broyden.accelerate())

dY = diff(rna.Xplus-rna.X,1,2);
dR = diff(rna.Xplus,1,2);

dim = size(dY,1);

%%

cvx_begin
    variable H(dim,dim)
    minimize norm(H*dR - dY,'fro') % + 1e-6*norm(H,'fro')
    subject to
    H == H'
%     H'*dR == dY
%     H*dR == dY
%     H >= 0
cvx_end

%%
reg = 1e-6*eye(size(dR,2));
H_test = -(dY/(dR'*dR+reg))*dR';

x = rna.X(:,end);
r = (rna.Xplus(:,end)-rna.X(:,end));
xplus = rna.Xplus(:,end);
x_bfgs = x+H*r;
x_test = x+H_test*r;
[finfo.f(x) finfo.f(rna.accelerate()) finfo.f(x_bfgs) finfo.f(x_test)]


