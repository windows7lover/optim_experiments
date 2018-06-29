% solve_svm

clearvars
close all
clc

% Algo param
nIte = 5000;
RmpeParam.lambda = 1e-1;
RmpeParam.lambdamin = 1e-15;


% svm param
npoints = 1000;
nfeatures = 200;
[X,y,xstar] = gen_problem_svm(npoints,nfeatures);
x0 = zeros(size(X,1),1);

%%
% Get Function
paramSvm.X = X;
paramSvm.y = y;
paramSvm.xstar = xstar;
paramSvm.x0 = x0;

finfo = getFunction('dualsvm',paramSvm);

% Regularization of finfo
epsilon = 1e-6;
D = norm(finfo.x0-finfo.xstar);
finforeg = finfo;
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
RmpeParam.k = 5;
RmpeParam.optialgo=@proximalgradient;
% RmpeParam.ampealgo=@ampe_box;
RmpeParam.ampealgo=@ampe;
RmpeParam.doLineSearch=true;
RmpeParam.useprox = true;
RmpeParam.doAdaptiveLambda = true;
RmpeParam.forceDecrease = true;
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

offset = 0*(1e-6);

figure
semilogy(ite_vec,abs(GradTolf(idx_algo)+offset),'-x','Color',colors(1,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
hold on
semilogy(ite_vec,abs(FistaTolf(idx_algo)+offset),'-*','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(ite_vec,abs(FistaRestartTolf(idx_algo))+offset,'-*','Color',colors(5,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(ite_vec,abs(RmpeTolf(idx_rmpe)+offset),'-o','Color',colors(3,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
% axis([-inf,inf,1e-14,inf])
legend({'gradient','Fista','Fista restarted','Acc Gradient','Acc Fista'})
ylabel('tolf','fontsize',16)
xlabel('iteration','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);



figure
semilogy(GradTime(idx_algo),abs(GradTolf(idx_algo)+offset),'-x','Color',colors(1,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
hold on
semilogy(FistaTime(idx_algo),abs(FistaTolf(idx_algo)+offset),'-*','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(FistaRestartTime(idx_algo),abs(FistaRestartTolf(idx_algo)+offset),'-*','Color',colors(5,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(RmpeTime(idx_rmpe),abs(RmpeTolf(idx_rmpe)+offset),'-o','Color',colors(3,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
% axis([-inf,inf,1e-14,inf])
legend({'gradient','Fista','Fista restarted','Acc Gradient','Acc Fista'})
ylabel('tolf','fontsize',16)
xlabel('iteration','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);


% save('lasso_problem1_solved')