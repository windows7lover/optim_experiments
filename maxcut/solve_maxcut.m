% Max cut

clearvars
close all
clc

nIte = 10000;
RmpeParam.lambda = 1e-6;
RmpeParam.lambdamin = 1e-20;


nNodes = 200;
nEdges = 2000;
G = randomGraph(nNodes,nEdges);
L = laplacian(G);
[xstar,ystar,fstar_primal,fstar_dual] = solve_maxcut_sdp(L);

param.laplacian = L;
param.xstar = ystar;
param.fstar= fstar_dual;
param.x0 = 0*ystar;

finfo = getFunction('maxcut',param);
%%

display('Subgradient')
param_subgrad.nite = 1;
[~,SubgradTolf,SubgradTolx,~,SubgradTime] = do_k_iterations(param_subgrad,@subgrad_method,nIte,finfo);
SubgradTime(end)



display('Dual averaging')
para_dual_ave.beta_hat = 1;
para_dual_ave.sk = 0*finfo.x0;
para_dual_ave.rho = sqrt(2*norm(finfo.x0-finfo.xstar)); % = sqrt(2D)
para_dual_ave.x0 = finfo.x0;
[~,DualaveTolf,DualaveTolx,~,DualaveTime] = do_k_iterations(para_dual_ave,@dual_averaging_method,nIte,finfo);
DualaveTime(end)

%%

display('RMPE Gradient')
RmpeParam.k = 25;
RmpeParam.optialgo=@subgrad_method;
RmpeParam.ampealgo=@ampe;
RmpeParam.doLineSearch=false;
RmpeParam.doAdaptiveLambda = false;
RmpeParam.forceDecrease = false;
nIteTemp = round(nIte/RmpeParam.k);
RmpeParam.optialgoparam=param_subgrad;
[~,SubgradRmpeTolf,SubgradRmpeTolx,~,SubgradRmpeTime] = do_k_iterations(RmpeParam,@abstract_ampe_adaptive_lambda_ls_stepsize,nIteTemp,finfo);
SubgradRmpeTime(end)

display('RMPE averaging')
RmpeParam.optialgo=@dual_averaging_method;
nIteTemp = round(nIte/RmpeParam.k);
RmpeParam.optialgoparam=para_dual_ave;
[~,DualaveRmpeTolf,DualaveRmpeTolx,~,DualaveRmpeTime] = do_k_iterations(RmpeParam,@abstract_ampe_adaptive_lambda_ls_stepsize,nIteTemp,finfo);
DualaveRmpeTime(end)

%% Plot


colors = [[150,150,150]; ...
[250,126,63]; ...
[115,123,13]; ...
[69,180,235]; ...
[235,111,217]];
colors = colors/255;

nPointsToPlot = 10;
idx_algo = round(linspace(1,length(SubgradTolf),nPointsToPlot+1)); idx_algo(1) = 1; idx_algo(end) = length(SubgradTolf);
idx_rmpe = round(linspace(1,length(SubgradRmpeTolf),nPointsToPlot+1)); idx_rmpe(1) = 1; idx_rmpe(end) = length(SubgradRmpeTolf);

offset = 0;

figure
semilogy(idx_algo,SubgradTolf(idx_algo),'-x','Color',colors(1,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
hold on
semilogy(1+(idx_rmpe-1)*RmpeParam.k,SubgradRmpeTolf(idx_rmpe)+offset,'-o','Color',colors(3,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
semilogy(idx_algo,DualaveTolf(idx_algo),'-x','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(1+(idx_rmpe-1)*RmpeParam.k,DualaveRmpeTolf(idx_rmpe)+offset,'-o','Color',colors(4,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
legend({'Subgradient','Acc subgradient','Dual averaging','Acc Dual averaging'})
ylabel('tolf','fontsize',16)
xlabel('iteration','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
axis tight



figure
semilogy(SubgradTime(idx_algo),SubgradTolf(idx_algo),'-x','Color',colors(1,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
hold on
semilogy(SubgradRmpeTime(idx_rmpe),SubgradRmpeTolf(idx_rmpe)+offset,'-o','Color',colors(3,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
semilogy(DualaveTime(idx_algo),DualaveTolf(idx_algo),'-x','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(DualaveRmpeTime(idx_rmpe),DualaveRmpeTolf(idx_rmpe)+offset,'-o','Color',colors(4,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
legend({'Subgradient','Acc subgradient','Dual averaging','Acc Dual averaging'})
ylabel('tolf','fontsize',16)
xlabel('iteration','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
axis tight