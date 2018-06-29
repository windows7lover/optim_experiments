

% savefig = true;
savefig = false;

if(savefig)
    dir = 'figs/';
    namefig_ite = [dir name '_nite'];
    namefig_time = [dir name '_time'];
end

%%

legendCell = {'Saga','SGD','SVRG','Katyusha','RMPE-Saga','RMPE-SGD','RMPE-SVRG','RMPE-Katyusha'};

colors = [[150,150,150]; ...
[250,126,63]; ...
[115,123,13]; ...
[69,180,235]; ...
[235,111,217]];
colors = colors/255;

%%

figure

semilogy(itersaga,pred_saga,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(1,:));
hold on;
semilogy(itersgd,pred_sgd,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(2,:));
semilogy(itersvrg,pred_svrg,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(3,:));
semilogy(iterkatyusha,pred_katyusha,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(4,:));
semilogy(iterrmpe_saga,pred_rmpe_saga,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(1,:));
semilogy(iterrmpe_sgd,pred_rmpe_sgd,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(2,:));
semilogy(iterrmpe_svrg,pred_rmpe_svrg,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(3,:));
semilogy(iterrmpe_katyusha,pred_rmpe_katyusha,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(4,:));
hold off;

xlabel('xlabel','fontsize',16)
ylabel('valf','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
legend(legendCell,'location','best')
axis tight


if(savefig)
    saveas(gcf,namefig_ite,'epsc')
    saveas(gcf,namefig_ite,'fig')
end


figure

semilogy(time_saga,pred_saga,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(1,:));
hold on;
semilogy(time_sgd,pred_sgd,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(2,:));
semilogy(time_svrg,pred_svrg,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(3,:));
semilogy(time_katyusha,pred_katyusha,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(4,:));
semilogy(time_rmpe_saga,pred_rmpe_saga,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(1,:));
semilogy(time_rmpe_sgd,pred_rmpe_sgd,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(2,:));
semilogy(time_rmpe_svrg,pred_rmpe_svrg,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(3,:));
semilogy(time_rmpe_katyusha,pred_rmpe_katyusha,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(4,:));
hold off;

xlabel('xlabel','fontsize',16)
ylabel('valf','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
legend(legendCell)
legend(legendCell,'location','best')
axis tight


if(savefig)
    saveas(gcf,namefig_time,'epsc')
    saveas(gcf,namefig_time,'fig')
end
