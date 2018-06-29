
% Note: In this file, we plot abs(f-fstar) to avoid negative values in the
% log. Usually, fstar has been computed using CVX up to precision 1e-8. In
% pratice, the precision is of the order ~1e-10.


%% General definitions (legend, colors)

legendCell = {'Grad','GradNest','GradNestBack',['RMPE',num2str(RmpeParam.k)],['RMPE',num2str(RmpeParam.k),'-LS']};

colors = [[150,150,150]; ...
[250,126,63]; ...
[115,123,13]; ...
[69,180,235]; ...
[235,111,217]];
colors = colors/255;


%% vectors for x axis
stepPlot1 = round(nIte/10);
stepPlot2 = round(stepPlot1/RmpeParam.k);

xGrad = [0:stepPlot1:nIte, nIte];
xrmpe5 = [0:stepPlot2:length(rmpeGradTolX5)-1,length(rmpeGradTolX5)-1];

time_grad = tgrad*xGrad/max(xGrad);
time_nest = tnest*xGrad/max(xGrad);
time_nestback = tnextback*xGrad/max(xGrad);
time_rmpe5 = trmpe5*xrmpe5/max(xrmpe5);
time_rmpe5_ls = trmpe5_ls*xrmpe5/max(xrmpe5);

%% Pot [k,f(x_k)]

figure


semilogy(xGrad,tolFgrad(1+xGrad),'-x','Color',colors(1,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
hold on;
semilogy(xGrad,abs(tolFnest(1+xGrad)),'-*','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(xGrad,abs(tolFnestback(1+xGrad)),'-*','Color',colors(5,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(5*xrmpe5,abs(rmpeGradTolF5(1+xrmpe5)),'-o','Color',colors(3,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
semilogy(5*xrmpe5,abs(rmpeGradTolF5_ls(1+xrmpe5)),'-s','Color',colors(4,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
hold off;

xlabel('Gradient Oracle Calls','fontsize',16)
ylabel('Duality Gap','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
legend(legendCell,'location','best')
axis tight
title(['Dataset: ', dataset, ', reg. param. ' , num2str(reg)],'fontsize',12)


figure

semilogy(time_grad,tolFgrad(1+xGrad),'-x','Color',colors(1,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
hold on;
semilogy(time_nest,abs(tolFnest(1+xGrad)),'-*','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(time_nestback,abs(tolFnestback(1+xGrad)),'-*','Color',colors(5,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');
semilogy(time_rmpe5,abs(rmpeGradTolF5(1+xrmpe5)),'-o','Color',colors(3,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
semilogy(time_rmpe5_ls,abs(rmpeGradTolF5_ls(1+xrmpe5)),'-s','Color',colors(4,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');
hold off;

xlabel('Time (sec)','fontsize',16)
ylabel('Duality Gap','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
legend(legendCell)
legend(legendCell,'location','best')
axis tight
title(['Dataset: ', dataset, ', reg. param. ' , num2str(reg)],'fontsize',12)
