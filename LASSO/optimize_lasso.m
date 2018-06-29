% solve_lasso

clearvars
% close all
clc

% Algo param
nIte = 100000;
RmpeParam.lambda = 1e-1;
RmpeParam.lambdamin = 1e-11;


% Lasso param
rho = 10; % magnitude of x*
lambda = 10; % Regularizer
d = 1000; % d>n
s = 50; % sparsity of the solution
n = 200;
cond_number = 1e2; % only for lasso_random

% Generate problem ||Ax-b||_2 + lambda ||x||_1
[A, b, xstar] = gen_problem_Lasso(rho, lambda, s, d, n);
% [A, b, xstar] = gen_problem_Lasso_random(n,d,s,cond_number,rho,lambda);
x0 = (A'*A + lambda*eye(d))\(A'*b);


% Get Function
paramLasso.A = A;
paramLasso.b = b;
paramLasso.xstar = xstar;
paramLasso.x0 = x0;
paramLasso.lambda = lambda;

finfo = getFunction('Lasso',paramLasso);
% tolfun.tolx = @(x) max(norm(x(finfo.xstar==0),1),1e-14); % norm 2
% tolfun.tolx = @(x) sum(x(finfo.xstar==0)~=0); % norm 0

% Regularization of finfo
epsilon = 1e-6;
D = norm(finfo.x0-finfo.xstar);
finforeg = finfo;
% finforeg.f = @(x) finfo.f(x) + epsilon*norm(x)^2/(2*D^2);
% finforeg.fp = @(x) finfo.fp(x) + epsilon*norm(x)/(D^2);

%%


display('Gradient')
[~,GradTolf,GradTolx,~,GradTime] = do_k_iterations([],@proximalgradient,nIte,finfo);
GradTime(end)
 
display('fista')
[~,FistaTolf,FistaTolx,~,FistaTime] = do_k_iterations([],@fista,nIte,finfo);
FistaTime(end)
 
display('fista with restart')
[~,FistaRestartTolf,FistaRestartTolx,~,FistaRestartTime] = do_k_iterations([],@fista_restart,nIte,finfo);
FistaRestartTime(end)

display('RMPE Gradient')
RmpeParam.k = 10;
RmpeParam.doLineSearch = true;
RmpeParam.doAdaptiveLambda = true;
RmpeParam.forceDecrease = true;
RmpeParam.optialgo=@proximalgradient;
RmpeParam.ampealgo=@ampe;
RmpeParam.useprox = false;
RmpeParam.stepsizeprox = 1/finfo.L;
nIteTemp = round(nIte/RmpeParam.k);
[~,RmpeTolf,RmpeTolx,~,RmpeTime] = do_k_iterations(RmpeParam,@abstract_ampe_adaptive_lambda_ls_stepsize,nIteTemp,finforeg);
RmpeTime(end)

% display('RMPE Fista')
% RmpeParam.optialgo=@fista;
% nIteTemp = round(nIte/RmpeParam.k);
% [~,RmpeFistaTolf,RmpeFistaTolx,~,RmpeFistaTime] = do_k_iterations(RmpeParam,@abstract_ampe_adaptive_lambda_ls_stepsize,nIteTemp,finforeg);
% RmpeFistaTime(end)

%% Plot


colors = [[150,150,150]; ...
[250,126,63]; ...
[115,123,13]; ...
[69,180,235]; ...
[235,111,217]];
colors = colors/255;

nPointsToPlot = 10;
ite_vec = linspace(0,nIte,nPointsToPlot+1);
idx_algo = 1+round(ite_vec);
idx_rmpe = 1+round(ite_vec/RmpeParam.k);

offset = 100*eps;

figure
semilogy(ite_vec,GradTolx(idx_algo),'-x','Color',colors(1,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
hold on
semilogy(ite_vec,FistaTolx(idx_algo)+offset,'-*','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(ite_vec,FistaRestartTolx(idx_algo)+offset,'-*','Color',colors(5,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(ite_vec,RmpeTolx(idx_rmpe)+offset,'-o','Color',colors(3,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
% semilogy(ite_vec,RmpeFistaTolf(idx_rmpe)+offset,'-s','Color',colors(4,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
% axis([-inf,inf,1e-14,inf])
legend({'gradient','Fista','Fista restarted','Acc Gradient','Acc Fista'})
ylabel('tolf','fontsize',16)
xlabel('iteration','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);


figure
semilogy(GradTime(idx_algo),GradTolx(idx_algo),'-x','Color',colors(1,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
hold on
semilogy(FistaTime(idx_algo),FistaTolx(idx_algo)+offset,'-*','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(FistaRestartTime(idx_algo),FistaRestartTolx(idx_algo)+offset,'-*','Color',colors(5,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(RmpeTime(idx_rmpe),RmpeTolx(idx_rmpe)+offset,'-o','Color',colors(3,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
% semilogy(RmpeFistaTime(idx_rmpe),RmpeFistaTolf(idx_rmpe)+offset,'-s','Color',colors(4,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
% axis([-inf,inf,1e-14,inf])
legend({'gradient','Fista','Fista restarted','Acc Gradient','Acc Fista'})
ylabel('tolf','fontsize',16)
xlabel('iteration','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);